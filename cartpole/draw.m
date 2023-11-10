function draw(y,t,params)

    % configuration coordiantes
    x = y(1);   % cart position
    th = y(2);  % pendulum angle
    
    M = params(1);
    m = params(2);
    L = params(3);

    % dimensions
    % L = 2;  % pendulum length
    W = 1*sqrt(M/5);  % cart width
    H = .5*sqrt(M/5); % cart height
    wr = .2; % wheel radius
    mr = .3*sqrt(m); % mass radius
    
    % positions
    % y = wr/2; % cart vertical position
    y = wr/2+H/2; % cart vertical position
    w1x = x-.9*W/2;
    w1y = 0;
    w2x = x+.9*W/2-wr;
    w2y = 0;
    
    % compute location of pole mass
    px = x + L*sin(th);
    py = y - L*cos(th);
    
    % plot the ground
    window_width = L+ 4*mr;
    plot([-500,500],[0,0],'w','LineWidth',2)
    hold on

    % draw cart
    % cart is official caltech orange
    rectangle('Position',[x-W/2,y-H/2,W,H],'Curvature',.1,'FaceColor',"#FF6C0C",'EdgeColor',[1 1 1])
    
    % draw wheels
    rectangle('Position',[w1x,w1y,wr,wr],'Curvature',1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])
    rectangle('Position',[w2x,w2y,wr,wr],'Curvature',1,'FaceColor',[1 1 1],'EdgeColor',[1 1 1])
    
    % draw pole
    plot([x px],[y py],'w','LineWidth',2)
    
    % draw pole mass
    rectangle('Position',[px-mr/2,py-mr/2,mr,mr],'Curvature',1,'FaceColor',"#595959",'EdgeColor',[1 1 1])

    % window properties
    xlim([x-window_width x+window_width]);
    axis equal
    set(gca,'Color','k','XColor','w','YColor','w')
    set(gcf,'Color','k')
    set(gcf,'InvertHardcopy','off')   
    
    % title
    str = sprintf("Time =  %.2f", t);
    title(str,'Color','w')

    box on
    drawnow
    set(gcf,'Position',[10 900 800 400])
    hold off
end
