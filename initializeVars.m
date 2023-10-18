function [A2,Cc,B1,B2,B3,T1,GNA,GNB,sA,sB,dsA,dsB,AeffA,AeffB,dt,Skip,NSteps,datestring,aA,aB,LeffA,LeffB] = initializeVars(A1,myRadA,myRadB,LA,LB,dt0,numOutput,totalTime)
    
    A2 = A1;  % bend modulus rotation 
    Cc = A1; % twist modulus about e3
    
    B1 = 40*A1; % Coordinate axis soft constraint
    B2 = 1*B1;
    B3 = 1*B1;
    T1 = 1*A1; % these params ensure coordinate axis does not rotate
    
    dsneededA = myRadA./3.5; % immersed boundary mod for radius (when using a 4 grid spread as originally shown by Peskin)
    dsneededB = myRadB./3.5; % immersed boundary mod for radius
    
    GNA = round(LA./dsneededA,0); % number of nodes defining cell A
    GNB = round(LB./dsneededB,0); % number of nodes defining cell B
    sA = linspace(0,LA,GNA); % first and last nodes are at the very end
    dsA = sA(2)-sA(1); % length between 2 nodes for cell "A"
    sB = linspace(0,LB,GNB); % first and last nodes are at the very end
    dsB = sB(2)-sB(1); % length between 2 nodes for cell "B"
    AeffA = dsA.*3.5;  % the factor of 3.5/2 (comes from the average of "k_w = 4,5,6, and another 4)
    AeffB = dsB.*3.5;  % the factor of 3.5/2 (comes from the average of "k_w = 4,5,6, and another 4)
    dt = dt0;
    samplet = totalTime./numOutput;    % how often to collect the data and plot
    Skip = (samplet)/dt; 
    NSteps = ceil(totalTime./dt./Skip); % total number of stored time steps
    datestring = datestr(now,'yy_mm_dd_HH_MM');
    aA = AeffA.*ones(GNA,1); % won't be used in code... only as a reference for plotting
    aB = AeffB.*ones(GNB,1); % won't be used in code... only as a reference for plotting
    LeffA = LA + 2*AeffA; % adjusted length
    LeffB = LB + 2*AeffB; % adjusted length
end