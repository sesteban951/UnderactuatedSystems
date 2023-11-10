%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CART-POLE SIM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all, close all, clc, clear figures;

% define system parameters
mc = 5;      % cart mass
mp = 2;      % pole mass
l = 2;       % pole length
g = 9.8;     % gravity
bc = 0;      % cart damping
bp = 0;      % pole damping

% pack into sys_info list
params(1) = mc;
params(2) = mp;
params(3) = l;
params(4) = g;
params(5) = bc;
params(6) = bp;

% evaluation and linearization functions
[f_x, g_x, Df, Dg] = make_funcs(params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATION

% simulation time
t_max = 8;

% initial 
x0 =  [-3;
       .9*pi; 
       -0;
       0]; % [pos, theta, vel, rot_vel] 

% final
xf =  [10;
       1*pi; 
       0;
       0]; % [pos, theta, vel, rot_vel] 

% x0 = manifold_x0(x0,params);

% MPC params
dt = 0.1;  % discretization size
N = 100;    % horizon length
nx = 4;    % state dimension size
nu = 1;    % input dimensions size

% total sizes for MPC
nx_tot = nx*(N-1);
nu_tot = nu*(N-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIMIZER
yalmip('clear')

% decision variables, super wide vectors as well
x_var = sdpvar(nx, N);      % states, [x1, ..., xN] a (4 X N) matrx
u_var = sdpvar(nu, N-1);    % inputs, [u1, ..., uN-1] (1xN-1)
  
x0_ = sdpvar(nx, 1);

Ad_ = sdpvar(nx, nx_tot);   % discrete drift
Bd_ = sdpvar(nx, nu_tot);               % discrete actuation
Cd_ = sdpvar(nx, N-1);                  % discrete residual

% containers to hold all linear discrete matrices, super wide matrices
Ad = zeros(nx, nx_tot);
Bd = zeros(nx, nu_tot);
Cd = zeros(nx, N-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRAINTS

% build constraint list
% state constraints, x in X
z_min = -200;
z_max =  200;
th_min = -100;
th_max = 100;

z_dot_min = -200;
z_dot_max =  200;
th_dot_min = -100;
th_dot_max = 100;

% input constraints, u in U
u_max =  1E3;
u_min = -1E4;

% matrices and vectors to describe polyhedron for state cons., Cx*x <= dx
% and input constraints, Cu*u <= du
Cx = kron(eye(nx),[1;-1]); % 8x4
dx = [z_max;     % put all constraints into polyhedron
     -z_min; 
     th_max; 
     -th_min; 
     z_dot_max; 
     -z_dot_min; 
     th_dot_max; 
     -th_dot_min];

Cu = [1;-1];
du = [u_max; -u_min];

% init constraint list and add the inital constraint to it
init_const = (x_var(:,1) == x0_);
constraints = [init_const];

% for loop to add all constraints
for i = 1 : N-1

    % new state constraint, described as polyhedron
    state_cons = (Cx*x_var(:,i) <= dx);

    % new input constraint
    input_cons = (Cu*u_var(:,i) <= du);

    % new dynamics constraint
    dyn_cons = x_var(:,i+1) == Ad_(:,(i-1)*nx+1 : i*nx)*x_var(:,i) ... 
               + Bd_(:,i)*u_var(i) ...
               + Cd_(:,i);
    
    % append constraints to the constraints list
    constraints = [constraints; input_cons; state_cons; dyn_cons];
end

% append terminal constraint
term_const = (x_var(:,N) == xf);
constraints = [constraints; term_const];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OBJECTIVE FUNCTION

% state and input scaling
x_scaling = [10, 50, 10, 10];
u_scaling = .01;

% build Qx and Qu matrices for the horizon in giant Q matrix
Qx = kron(eye(N), diag(x_scaling));
Qu = kron(eye(N-1), diag(u_scaling));
Q = blkdiag(Qx, Qu);  % [Qx1, ..., QxN, Qu1, ..., QuN-1]

% terminal cost
term_scaling = [1,1,1,1];
V = diag(term_scaling);

% build objective, stack x, u -> J = [x; u]' * Q * [x;u] 
x_var_e = x_var(:) - kron(ones(N,1),xf);
objective = (1/2)*([x_var_e; u_var(:)]' * Q * [x_var_e; u_var(:)]...
            + (x_var(:,N)-xf)'*V*(x_var(:,N)-xf));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIMIZER

% configure optimizer settings
options = sdpsettings('solver','mosek','verbose',0);

% create the optimzer object/function
P = optimizer(constraints, ...                 % constraints
              objective, ...                   % objective
              options,...                      % options
              {Ad_, Bd_, Cd_, x0_},...         % inpute paramters
              {x_var,u_var});                % dec. vars.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use LQR to compute fixed control policy, u = -K_lqr * x
% linearize about the desired target, x_dot = A_xf x + B_xf u
A_xf = Df(xf(1),xf(2),xf(3),xf(4));
B_xf = g_x(xf(1),xf(2),xf(3),xf(4));

% cost matrices
Qx_lqr = diag(x_scaling);
Qu_lqr = diag(u_scaling);

% will use this gain matrix for forward integarting 
K_lqr = lqr(A_xf, B_xf, Qx_lqr, Qu_lqr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% containers for storing stuff
T = [0];    % history of time steps
X = [x0'];  % history of states at time steps
U = [];     % history of inputs

% contaniers to store MPC trajectories
X_MPC_Path = zeros (N, nx, ceil(t_max/dt));
U_MPC_Path = zeros (N-1, nu, ceil(t_max/dt));
T_MPC_Path = zeros(N,1,ceil(t_max/dt));
x_path_ind = 1;

% fwd simulate for every time step
wbar = waitbar(0,'Simulating');
t_now=0;
for t = 0:dt:t_max
    
    waitbar((t_now/t_max),wbar, "LQR Fwd Sim for ref traj.")
    % fwd simulate for small dt
    t_lqr = 0:dt:dt*(N-1);  % horizon length
    [~,x_] = ode45(@(t,x) dynamics_fbk(t,x,params,f_x,g_x, K_lqr,xf), ... 
                         t_lqr, x0);
    
    % reulting traj w/ LQR fwd sim
    x_bar = x_' ;
    x_eq = x_'-kron(ones(1,N),xf);
    u_bar = -K_lqr*x_eq; %  computing u_k = -K * x_k

    % Build Ad, Bd, Cd for every integration step
    Ad = [];
    Bd = [];
    Cd = [];

    waitbar((t_now/t_max),wbar, "Building Ak, Bk, Ck.")
    % linearize about the reference trajectories
    for k = 1:N-1
        % which x_bar and u_bar
        x_k = x_bar(:,k);
        u_k = u_bar(:,k);

        % compute f(x), g(x) linearizations
        Df_k = Df(x_k(1),x_k(2), x_k(3), x_k(4)); 
        Dg_k = Dg(x_k(1),x_k(2), x_k(3), x_k(4));  
        f_x_k = f_x(x_k(1),x_k(2), x_k(3), x_k(4));
        g_x_k = g_x(x_k(1),x_k(2), x_k(3), x_k(4));

        % continous time linear matrices
        Ac = Df_k + Dg_k * u_k;
        Bc = g_x_k;   
        Cc = f_x_k + g_x_k * u_k - Ac * x_k - Bc * u_k;

        % continous to discrete time linear matrices
        [Ad_k,Bd_k,Cd_k] = css2dss('Exact',dt,Ac,Bc,Cc);  %%%%%%%%%%%%%

        % append to discrete matrix list
        Ad = [Ad , Ad_k];
        Bd = [Bd , Bd_k];
        Cd = [Cd , Cd_k];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this is where the magic happens -- optimzer plus dynamics solver
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % solve MPC problem for x_star and u_star given init state x_0
    % solution will return optimal trajectory and optimal input
    waitbar((t_now/t_max),wbar, "Solving OCP.")
    [sol, diagnostics,d1,d2,d3,d4] = P({Ad,Bd,Cd,x0});
    x_star = sol{1};    
    u_star = sol{2};

    if (diagnostics ~= 0)
        disp("Not Feasible at time:")
        disp(t_now)
        close(wbar)
        break
    end

    % take first input and use as feed forward via zero order hold
    x_des = x_star(:,1);
    u_ff = u_star(:,1);
    % Note: you can also devise some kind of low level control to track
    % the optimal state, x_star. Can do it here or inside ODE dynamics

    % fwd simulate with given x0 and u_ff to see where you end up after the
    % small time step

    % [0,dt] is the span of time in which you want ot integrate
    % x_0 is the intial condition for the solver
    % PendulumODE is the dynamics of your system with inputs 
    [t,x] = ode45(@(t,x) dynamics_ff(t,x,params,f_x,g_x,x_des,u_ff), [0 dt], x0);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % store time steps and states for actual evolution of system
    t_now = T(end);
    

    T = [T; t+T(end)]; % add new batch of time steps, add T(end) to
                       % splice and make continous stream of time
    X = [X; x];        % add new batch of states at time steps
    U = [U; ones(length(t), 1) * u_ff]; % add new batch of inputs

    % store MPC optimal inputs and states
    x_star = x_star';

    % store MPC trajectories at each step (Noel did this part)
    X_MPC_Path(:,:,x_path_ind) = x_star;
    U_MPC_Path(:,:,x_path_ind) = u_star;
    T_MPC_Path(:,:,x_path_ind) = (t_now+(0:dt:dt*(N-1)))';
    x_path_ind = x_path_ind+1;

    % new inital condition is the last state after MPC u_ff is applied
    x0 = x(end,:)';

end
waitbar((t_now/t_max),wbar, "Finished Sim.")
close(wbar)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PLOT
% plot state information

labels = ["Cart Pos., $z$";
          "Pole Pos., $\theta$";
          "Cart Vel., $\dot{z}$";
          "Pole Vel., $\dot{\theta}$"];
% 
% figure(1);
% for i = 1:length(x0)
%     subplot(2,2,i)
%     plot(T,X(:,i))
%     set(gca,'Color','w','XColor','w','YColor','w')
%     set(gcf,'Color','k')
%     yline(xf(i))
%     xlabel("Time, $t$",'Interpreter','latex','Color','white')
%     ylabel(labels(i),'Interpreter','latex','Color','white')
%     grid on
% end

% draw the cart pole simulation
%%
figure(2);
size_des = 600;     % has to be less than total data size
n = size(T)/size_des;
n = round(n);
T_ = T(1:n:end);
X_ = X(1:n:end,:);
drawnow
tic
while i <= length(T_)
    while toc < T_(i)
        %
    end
    draw(X_(i,:),T_(i), params)
    i = i+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AUXILLARY FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dynamics
function x_dot = dynamics_ff(~,x,params,f_x,g_x,x_des,u_ff)
    
    % PD controller
    u_PD = - [1,0]*[x(1)-x_des(1); x(3)-x_des(3)] -...
            [1,0]*[x(2)-x_des(2); x(4)-x_des(4)];
    u_PD = 0;
    u = u_PD + u_ff;

    % set dynamics info
    f = f_x(x(1),x(2),x(3),x(4));
    g = g_x(x(1),x(2),x(3),x(4));
    x_dot = f+ g*u;
end

% dynamics with LQR for fwd prop
function x_dot = dynamics_fbk(~,x,params,f_x,g_x, K_lqr,xf)

    % control input
    u = -K_lqr*(x-xf);

    % set dynamics info
    f = f_x(x(1),x(2),x(3),x(4));
    g = g_x(x(1),x(2),x(3),x(4));
    x_dot = f+ g*u;

end


% get linearizations of vectors, Df, Dg at point x_k
function [f_x,g_x,Df, Dg] = make_funcs(params)

    syms th th_dot z z_dot

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % extract system parameters
    mc = params(1);  % cart mass
    mp = params(2);  % pole mass
    l = params(3);  % pole length
    g = params(4);  % gravity
    bc = params(5); % cart damping
    bp = params(6); % pole damping

    % convert to configuration space
    q_dot = [z_dot; th_dot];

    % Matrix info
    D = [(mc+mp),      mp*l*cos(th); 
         mp*l*cos(th), mp*(l^2)];
    
    C = [bc, -mp*l*th_dot*sin(th); 
          0, bp];
    
    G = [0; 
         mp*g*l*sin(th)];
    
    B = [1; 
         0];
    
    % set dynamics info
    f_x = [q_dot;
         -inv(D)*(C*q_dot + G)];
    g_x = [zeros(2,1);
          inv(D)*B]; 

    % compute jacobians
    Df = [diff(f_x,z), diff(f_x,th), diff(f_x,z_dot), diff(f_x,th_dot)];
    Dg = [diff(g_x,z), diff(g_x,th), diff(g_x,z_dot), diff(g_x,th_dot)];

    % state values to compute
    f_x = matlabFunction(f_x, 'Vars',{'z','th','z_dot','th_dot'});
    g_x = matlabFunction(g_x, 'Vars',{'z','th','z_dot','th_dot'});
    Df = matlabFunction(Df, 'Vars',{'z','th','z_dot','th_dot'});
    Dg = matlabFunction(Dg, 'Vars',{'z','th','z_dot','th_dot'});

end

% compute random crap after you're done
function [L,Pz,x_] = compute_stuff(x,params)

    % extract system parameters
    mc = params(1);  % cart mass
    mp = params(2);  % pole mass
    l = params(3);  % pole length
    g = params(4);  % gravity
    bc = params(5); % cart damping
    bp = params(6); % pole damping

    [row,col] = size(x);

    % extract state information
    z  = x(:,1);
    th = x(:,2);
    z_dot = x(:,3);
    th_dot= x(:,4);
    
    Pz = [];
    L = [];
    for i = 1:row
        % compute momentum, this thing should be constatn for all time
        % assuming no damping?
        Pz_ = (mc+mp)*z_dot(i) + mp*th_dot(i)*l*cos(th(i));
        Pz = [Pz Pz_];

        % compute Lagrangian
        L_ = 0.5*(mc+mp)*(z_dot(i)^2) + mp*z_dot(i)*th_dot(i)*l*cos(th(i)) + ...
            0.5*mp*(l^2)*(th_dot(i)^2) + mp*g*l*cos(th(i));
        L = [L L_];
    end
    
    x_ = x;

end

