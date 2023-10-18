%% Unconstrained Filament Algorithm... Based on Peskin,Lim,& SD Olson's work March 7, 2022
function [tInfo] = StackedCellsFluidForce_IB(LA,LB,myRadA,myRadB)

%% Algorithm Notes
% Code assumes a flat wall (floor) with two cylindrical on top of it. 
% These two cells are stacked on top of eachother and the fluid velocity is calculated using the immersed boundary method (Peskin)
% Currently the cells are fixed in place in their initialized configuration

% The force from the fluid velocity on the top cell can be extracted from
% the data to establish a baseline for rupture strength in bacterial
% adhesin experiments

% Uses a 2nd order Runge Kutta (RK) 
% Uses Cosserat model for coordinate axis (soft constraint) to mimic sheer on the cell wall
% initializes the fluid flow everywhere in Y then allows the forces from
% cells to equilibrate and create true stokes flow pattern (Easy update to include turbulence)
%% User Input Variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
myVersion = 'Version1_StraightCellsOnly';

%% Cell inputs 
LA = 2.0;       % length in microns. Cell "A" (effective length will have 2*radius added to this number)
LB = 2.0;       % length in microns. Cell "B"
myRadA = 0.5;   % radius of each E.Coli cell (smaller radii need smaller dt)
myRadB = 0.5;   % radius of each E.Coli cell (smaller radii need smaller dt)

yshift = 0.0;   % y shift from center of bottom TO center of top [um]
zshift = -0.4;  % z shift from center of bottom TO center of top [um]
phiA = pi./6;   % angle of BactA from z hat
phiB = -pi./8;  % angle of BactB from z hat

kappa1 = 0.0;   % These MUST be zero for now... may be modified to non-zero if INITIAL shape and location generation are modified
kappa2 = 0.0;
tau0 = 0.0; 
A1 = 10.0;      % bend modulus about e1 (assumes each cell has the same stiffness and are symmetric in e1 and e2)
%% Time inputs
totalTime = 0.02;   % total time the program will run [seconds] 
dt0 = 4E-5;         % smallest time step. 4E-5 work [seconds]
numOutput = 200;    % number of times the cells and velocity is plotted during convergnce

%% Environment inputs
uInit = 10.0;       % velocity of fluid in y hat direction in [um/s] (integrating fluid per volume at boundary will give true fluid velocity / time)
eta = 0.001;        % viscosity [Pa*s] (0.001 = water viscosity)
rho = 1E-6;         % density of fluid in [ug/(um^3)] (water = 1E-6 in our units)

%% plotting inputs
colorvecta = [0.1 0.8 0.8];     % [RGB] color of bact A
colorvectb = [0.8 0.1 0.2];     % [RGB] color of bact B

% establishes the angle at which to view the two generated plots
mySideView1 = 90; 
mySideView2 = -85;
myTopView1 = 90;
myTopView2 = -30;
% establishes which eulerian grid points will be plotted [1:skip:end]
mySkipX = 2;
mySkipYZ = 2;

% end of user inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize More Cell and Timing Variables
[A2,Cc,B1,B2,B3,T1,GNA,GNB,sA,sB,dsA,dsB,AeffA,AeffB,dt,Skip,NSteps,datestring,aA,aB,LeffA,LeffB] = initializeVars(A1,myRadA,myRadB,LA,LB,dt0,numOutput,totalTime);

%% Immersed Boundary Fluid Initialization
[boxSize,hFL,fLNum,ipIB,imIB,ygrid,xgrid,zgrid] = initFluid(LA,LB,dsA,dsB);

% set up fluid flow
uFL = zeros(fLNum,fLNum,fLNum,3); % Initial fluid velocity 
uFL(:,1,:,2) = uInit;
%%
aDelta = init_aDelta_3D(fLNum,dt,eta,hFL,rho); % aDelta variable for the IB delta function

%% Timing
Time = zeros(NSteps,1); 

%% Define Derivatives 6th order if there's enough points along cell

if GNA>6
    L1A = FirstDer6(GNA,dsA);
else
    L1A = sparse((2:GNA-1),(1:GNA-2),-0.5,GNA,GNA) ...
   + sparse((2:GNA-1),(3:GNA),0.5,GNA,GNA) ...
   + sparse((1),([1 2]),[-1 1],GNA,GNA) ...
   + sparse((GNA),([GNA-1 GNA]),[-1 1],GNA,GNA);
end
if GNB>6
    L1B = FirstDer6(GNB,dsB);
else
    L1B = sparse((2:GNB-1),(1:GNB-2),-0.5,GNB,GNB) ...
   + sparse((2:GNB-1),(3:GNB),0.5,GNB,GNB) ...
   + sparse((1),([1 2]),[-1 1],GNB,GNB) ...
   + sparse((GNB),([GNB-1 GNB]),[-1 1],GNB,GNB);
end

% averaging matrices
avMatA = sparse((1:GNA-1),(1:GNA-1),0.5,GNA,GNA-1) ...
   + sparse((2:GNA),(1:GNA-1),0.5,GNA,GNA-1);
avMatB = sparse((1:GNB-1),(1:GNB-1),0.5,GNB,GNB-1) ...
   + sparse((2:GNB),(1:GNB-1),0.5,GNB,GNB-1);

%% Create relaxed bend and twist moduli
k1A = kappa1.*ones(GNA,1); % can't make these non-zero for now
k2A = kappa2.*ones(GNA,1);
tauA = tau0.*ones(GNA,1);

k1B = kappa1.*ones(GNB,1); % can't make these non-zero for now
k2B = kappa2.*ones(GNB,1);
tauB = tau0.*ones(GNB,1);
%% Create cell locations and orthonormal axes
[ra,rb,xa,ya,za,xb,yb,zb,d1a,d1b,d2a,d2b,d3a,d3b,rOriga,rOrigb] = initLocs(aA,aB,sA,sB,phiA,phiB,GNA,GNB,boxSize,AeffA,AeffB,yshift,zshift);
%% Create wall location and orthonormal axes
[AeffW,wallR,LmxMin,LmxMax,LmyMin,LmyMax,LmzMin,LmzMax,d1W,d2W,d3W,cntW] = createWall(fLNum,dsA,dsB,ya,za,yb,zb,boxSize,xa,xb,AeffA,AeffB);
wallROrig = wallR; % sets up code so that wall doesn't move

%% A matrix
[d1Halfa,d2Halfa,d3Halfa] = AMatHalfRot(d1a,d2a,d3a,GNA); % half angle to start us off
[d1Halfb,d2Halfb,d3Halfb] = AMatHalfRot(d1b,d2b,d3b,GNB);

%% video title
vidtitle1 = strcat('Peptide_',myVersion,'_',datestring,'_eta',num2str(eta*1000, '%6.3g'),'_GNA',num2str(GNA, '%6.3g'),'_dt', num2str(dt, '%6.3g'),'_A',num2str(A1,'%6.3g'),'_B',num2str(B1,'%6.3g'),'_L',num2str(LeffA));
vidtitle2 = strcat('PeptideSide_',myVersion,'_',datestring,'_eta',num2str(eta*1000, '%6.3g'),'_GNA',num2str(GNA, '%6.3g'),'_dt', num2str(dt, '%6.3g'),'_A',num2str(A1,'%6.3g'),'_B',num2str(B1,'%6.3g'),'_L',num2str(LeffA));

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
for jj=2:NSteps % will plot on each jj = X, k=1
 
    xa(:,jj) = xa(:,jj-1);ya(:,jj) = ya(:,jj-1);za(:,jj) = za(:,jj-1);
    xb(:,jj) = xb(:,jj-1);yb(:,jj) = yb(:,jj-1);zb(:,jj) = zb(:,jj-1);
    Time(jj) = Time(jj-1);

for k = 1:Skip
for rk = 1:2      %% Runge-Kutta loop
    uFL(:,1,:,2) = uInit;
    dt = dt0./2;
    Time(jj) = Time(jj) + dt;

%% All lagrangian calculations and updates in here

k1halfa = (k1A(2:GNA) + k1A(1:GNA-1))./2;
k2halfa = (k2A(2:GNA) + k2A(1:GNA-1))./2;
tauhalfa = (tauA(2:GNA) + tauA(1:GNA-1))./2;

k1halfb = (k1B(2:GNB) + k1B(1:GNB-1))./2;
k2halfb = (k2B(2:GNB) + k2B(1:GNB-1))./2;
tauhalfb = (tauB(2:GNB) + tauB(1:GNB-1))./2;

% defines the first derivative of node locations using central difference in Olson paper
[xsCDa,ysCDa,zsCDa] = FirstDerCD1d(xa(:,jj),ya(:,jj),za(:,jj),dsA); 
[xsCDb,ysCDb,zsCDb] = FirstDerCD1d(xb(:,jj),yb(:,jj),zb(:,jj),dsB); 

% defines the first derivative of d basis using central difference in Olson paper (3 dimensions are x,y,z)
[d1sCDa,d2sCDa,d3sCDa] = FirstDerCD3d(d1a,d2a,d3a,dsA);
[d1sCDb,d2sCDb,d3sCDb] = FirstDerCD3d(d1b,d2b,d3b,dsB);


[d2sCDDotd3Halfa] = dotUnit(d2sCDa,d3Halfa); % N1half term for Eq.24
[d3sCDDotd1Halfa] = dotUnit(d3sCDa,d1Halfa); % N2half term
[d1sCDDotd2Halfa] = dotUnit(d1sCDa,d2Halfa); % N3half term

[d2sCDDotd3Halfb] = dotUnit(d2sCDb,d3Halfb); % N1half term for Eq.24
[d3sCDDotd1Halfb] = dotUnit(d3sCDb,d1Halfb); % N2half term
[d1sCDDotd2Halfb] = dotUnit(d1sCDb,d2Halfb); % N3half term
    
d1FHalfa = d1Halfa(:,1).*xsCDa + d1Halfa(:,2).*ysCDa + d1Halfa(:,3).*zsCDa;
d2FHalfa = d2Halfa(:,1).*xsCDa + d2Halfa(:,2).*ysCDa + d2Halfa(:,3).*zsCDa;   
d3FHalfa = d3Halfa(:,1).*xsCDa + d3Halfa(:,2).*ysCDa + d3Halfa(:,3).*zsCDa;

d1FHalfb = d1Halfb(:,1).*xsCDb + d1Halfb(:,2).*ysCDb + d1Halfb(:,3).*zsCDb;
d2FHalfb = d2Halfb(:,1).*xsCDb + d2Halfb(:,2).*ysCDb + d2Halfb(:,3).*zsCDb;   
d3FHalfb = d3Halfb(:,1).*xsCDb + d3Halfb(:,2).*ysCDb + d3Halfb(:,3).*zsCDb;
%% Force at the half nodes
F1Halfa = B1.*(d1FHalfa); % scalar in the d1_half (k+1/2) direction
F2Halfa = B2.*(d2FHalfa);
F3Halfa = B3.*(d3FHalfa-1);

F1Halfb = B1.*(d1FHalfb); % scalar in the d1_half (k+1/2) direction
F2Halfb = B2.*(d2FHalfb);
F3Halfb = B3.*(d3FHalfb-1.0);

%% Moment of each bacterium
N1Halfa = A1.*(d2sCDDotd3Halfa - k1halfa);  % in the d1_half direction
N2Halfa = A2.*(d3sCDDotd1Halfa - k2halfa);
N3Halfa = Cc.*(d1sCDDotd2Halfa - tauhalfa);

N1Halfb = A1.*(d2sCDDotd3Halfb - k1halfb);  % in the d1_half direction
N2Halfb = A2.*(d3sCDDotd1Halfb - k2halfb);
N3Halfb = Cc.*(d1sCDDotd2Halfb - tauhalfb);

    %% Then assemble the half forces and moments as vectors IN CARTESIAN SPACE
    FVectHalfa = F1Halfa.*d1Halfa + F2Halfa.*d2Halfa + F3Halfa.*d3Halfa;
    NVectHalfa = N1Halfa.*d1Halfa + N2Halfa.*d2Halfa + N3Halfa.*d3Halfa;

    FVectHalfb = F1Halfb.*d1Halfb + F2Halfb.*d2Halfb + F3Halfb.*d3Halfb; % in cartesian space now
    NVectHalfb = N1Halfb.*d1Halfb + N2Halfb.*d2Halfb + N3Halfb.*d3Halfb;

    %% Then find force/length @k and moment/length @k as [F_{k+1/2}-F_{k-1/2}]/ds pg 284 of Lim et al Closed rod dynamics SIAM 2008
    % and similarly for moment per length
    fpLA = zeros(GNA,3);
    fpLA(1,:) = FVectHalfa(1,:)./dsA;  
    fpLA(GNA,:) = -FVectHalfa(end,:)./dsA;  
    [fpLmidA,~,~] = FirstDerCD3d(FVectHalfa,0,0,dsA);% defines the first derivative of d basis using central difference in Olson paper
    fpLA(2:GNA-1,:) = fpLmidA; % equation 30 in Olson paper... this is force per unit length at kth node

    fpLB = zeros(GNB,3);
    fpLB(1,:) = FVectHalfb(1,:)./dsB;  
    fpLB(GNB,:) = -FVectHalfb(end,:)./dsB;  
    [fpLmidB,~,~] = FirstDerCD3d(FVectHalfb,0,0,dsB);% defines the first derivative of d basis using central difference in Olson paper
    fpLB(2:GNB-1,:) = fpLmidB; % equation 30 in Olson paper... this is force per unit length at kth node


    r3Da = [xa(:,jj),ya(:,jj),za(:,jj)];
    r3Db = [xb(:,jj),yb(:,jj),zb(:,jj)];

    %% Writing out first term of equation 31 separate first... then cross product
    % boundary conditions are: force and moment at 1/2 and GridNum + 1/2 are all zero
    npLA = zeros(GNA,3);
    npLA(1,:) = NVectHalfa(1,:)./dsA;  
    npLA(GNA,:) = -NVectHalfa(end,:)./dsA;  
    [npLmidA,~,~] = FirstDerCD3d(NVectHalfa,0,0,dsA);% defines the first derivative of d basis using central difference in Olson paper
    npLA(2:GNA-1,:) = npLmidA; % equation 31 in Olson paper... this is moment per unit length at kth node (only have N terms)

    npLB = zeros(GNB,3);
    npLB(1,:) = NVectHalfb(1,:)./dsB;  
    npLB(GNB,:) = -NVectHalfb(end,:)./dsB;  
    [npLmidB,~,~] = FirstDerCD3d(NVectHalfb,0,0,dsB);% defines the first derivative of d basis using central difference in Olson paper
    npLB(2:GNB-1,:) = npLmidB; % equation 31 in Olson paper... this is moment per unit length at kth node (only have N terms)


    %% New way to find cross product npL terms 
    r3DCDa = [xsCDa,ysCDa,zsCDa];
    [npLNewA(:,1),npLNewA(:,2),npLNewA(:,3)] = crossProduct(r3DCDa,FVectHalfa);
    npLXTotalA = [avMatA*npLNewA(:,1),avMatA*npLNewA(:,2),avMatA*npLNewA(:,3)];
    npLA = npLA + npLXTotalA;

    r3DCDb = [xsCDb,ysCDb,zsCDb];
    [npLNewB(:,1),npLNewB(:,2),npLNewB(:,3)] = crossProduct(r3DCDb,FVectHalfb);
    npLXTotalB = [avMatB*npLNewB(:,1),avMatB*npLNewB(:,2),avMatB*npLNewB(:,3)];
    npLB = npLB + npLXTotalB;

%% spread the force on the fluid. This is ALL in cartesian coords
ra = [xa(:,jj),ya(:,jj),za(:,jj)];
rb = [xb(:,jj),yb(:,jj),zb(:,jj)];

% Bact A
FatNodeA = fpLA.*dsA;
FatNodeA = FatNodeA + T1.*(rOriga-ra);
MatNodeA = npLA.*dsA;
MRote3A = T1.*(tan(d1a(:,2)./d1a(:,1))).*d3a; % this term ensures bacteria A (Bottom) doesn't rotate
MatNodeA = MatNodeA - MRote3A;

% Bact B
FatNodeB = fpLB.*dsB;
FatNodeB = FatNodeB + T1.*(rOrigb-rb);
MatNodeB = npLB.*dsB;
MRote3B = T1.*(tan(d1b(:,2)./d1b(:,1))).*d3b; % this term ensures bacteria B (top) doesn't rotate
MatNodeB = MatNodeB - MRote3B;   

% Wall
FatNodeW = T1.*(wallROrig-wallR);
MRote3W = T1.*(tan(d1W(:,2)./d1W(:,1))).*d3W; % this term ensures bacteria B (top) doesn't rotate
MatNodeW =  -MRote3W;   

if rk == 1
    uFLatT = uFL;
    dt = dt0./2;
    
    wallRatT = wallR;
    xatTa = xa(:,jj); yatTa = ya(:,jj); zatTa = za(:,jj);
    d1atTa = d1a; d2atTa = d2a; d3atTa = d3a;
    xatTb = xb(:,jj); yatTb = yb(:,jj); zatTb = zb(:,jj);
    d1atTb = d1b; d2atTb = d2b; d3atTb = d3b;

    M4DperVA = vec_spread3D(MatNodeA,r3Da,hFL,fLNum,GNA); % Moment / volume
    M4DperVB = vec_spread3D(MatNodeB,r3Db,hFL,fLNum,GNB); % Moment / volume
    M4DperVW = vec_spread3D(MatNodeW,wallR,hFL,fLNum,cntW); % Moment / volume % maybe should use r3datT
    M4DperV = M4DperVA+M4DperVB+M4DperVW;
    
    g0XMompV = g0crossFun(M4DperV,hFL,ipIB,imIB); % Force / volume
    fpuvOnFluidA = vec_spread3D(FatNodeA,r3Da,hFL,fLNum,GNA); % gives a force per unit volume on the fluid of density rho and viscosity mu (mu=eta)
    fpuvOnFluidB = vec_spread3D(FatNodeB,r3Db,hFL,fLNum,GNB); % gives a force per unit volume on the fluid of density rho and viscosity mu (mu=eta)
    fpuvOnFluidW = vec_spread3D(FatNodeW,wallR,hFL,fLNum,cntW); % gives a force per unit volume on the fluid of density rho and viscosity mu (mu=eta)
   
    fpuvOnFluid = fpuvOnFluidA + fpuvOnFluidB + fpuvOnFluidW + g0XMompV;

    [uFL] = fluid3D_IB(uFL,fpuvOnFluid,aDelta,dt,eta,ipIB,imIB,hFL,rho); % Step Fluid Velocity
    
    uFLinterpA = vec_interp3D(uFL,r3Da,hFL,GNA,fLNum); % interpolated force at my nodes. 
    uFLinterpB = vec_interp3D(uFL,r3Db,hFL,GNB,fLNum); % interpolated force at my nodes. 
    uFLinterpW = vec_interp3D(uFL,wallR,hFL,cntW,fLNum); % interpolated force at my nodes. 
    
    g0XuFL = g0crossFun(uFL,hFL,ipIB,imIB);
    
    omegaA = vec_interp3D(g0XuFL,r3Da,hFL,GNA,fLNum); % should give an x,y,z component vector of torque
    omegaB = vec_interp3D(g0XuFL,r3Db,hFL,GNB,fLNum); % should give an x,y,z component vector of torque
    omegaW = vec_interp3D(g0XuFL,wallR,hFL,cntW,fLNum); % should give an x,y,z component vector of torque

    %%   
    omegaMagA = sqrt(omegaA(:,1).^2 + omegaA(:,2).^2 + omegaA(:,3).^2);
    omegaScA = omegaA./omegaMagA;
    omegaMagB = sqrt(omegaB(:,1).^2 + omegaB(:,2).^2 + omegaB(:,3).^2);
    omegaScB = omegaB./omegaMagB;
    omegaMagW = sqrt(omegaW(:,1).^2 + omegaW(:,2).^2 + omegaW(:,3).^2);
    omegaScW = omegaW./omegaMagW;
    
    %% Update positions (1/2 step forward)
        xa(:,jj) = xa(:,jj) + uFLinterpA(:,1).*dt;
        ya(:,jj) = ya(:,jj) + uFLinterpA(:,2).*dt;
        za(:,jj) = za(:,jj) + uFLinterpA(:,3).*dt;

        xb(:,jj) = xb(:,jj) + uFLinterpB(:,1).*dt;
        yb(:,jj) = yb(:,jj) + uFLinterpB(:,2).*dt;
        zb(:,jj) = zb(:,jj) + uFLinterpB(:,3).*dt;

        wallR = wallR + uFLinterpW.*dt;    
       
        %% Update the d hat bases (1/2 step)
    [d1newA,d2newA,d3newA] = updateDHat(GNA,dt,omegaMagA,d1a,d2a,d3a,omegaScA); % cell A
    [d1newB,d2newB,d3newB] = updateDHat(GNB,dt,omegaMagB,d1b,d2b,d3b,omegaScB); % cell B
    [d1newW,d2newW,d3newW] = updateDHat(cntW,dt,omegaMagW,d1W,d2W,d3W,omegaScW); % wall

     d1a = d1newA;d2a = d2newA;d3a = d3newA;
     d1b = d1newB;d2b = d2newB;d3b = d3newB;
     d1W = d1newW;d2W = d2newW;d3W = d3newW;
     
     %%
    [d1Halfa,d2Halfa,d3Halfa] = AMatHalfRot(d1a,d2a,d3a,GNA);
    [d1Halfb,d2Halfb,d3Halfb] = AMatHalfRot(d1b,d2b,d3b,GNB);

else % when rk == 2 (this is calculating full step with all t+dt/2 variables and applies everything to time = t pointss, fluid, rotations, etc)
    dt = dt0;
    r3Da = [xa(:,jj),ya(:,jj),za(:,jj)];
    r3DatTA = [xatTa,yatTa,zatTa];

    r3Db = [xb(:,jj),yb(:,jj),zb(:,jj)];
    r3DatTB = [xatTb,yatTb,zatTb];

    %% Doing this all at r3DatT
    M4DperVA = vec_spread3D(MatNodeA,r3DatTA,hFL,fLNum,GNA); % Moment / volume % maybe should use r3datT
    M4DperVB = vec_spread3D(MatNodeB,r3DatTB,hFL,fLNum,GNB); % Moment / volume % maybe should use r3datT
    M4DperVW = vec_spread3D(MatNodeW,wallRatT,hFL,fLNum,cntW); % Moment / volume % maybe should use r3datT
    M4DperV = M4DperVA + M4DperVB + M4DperVW;

    g0XMompV = g0crossFun(M4DperV,hFL,ipIB,imIB); % Force / volume because of moment
    
    fpuvOnFluidA = vec_spread3D(FatNodeA,r3DatTA,hFL,fLNum,GNA); % gives a force per unit volume on the fluid of density rho and viscosity mu (mu=eta)
    fpuvOnFluidB = vec_spread3D(FatNodeB,r3DatTB,hFL,fLNum,GNB); % gives a force per unit volume on the fluid of density rho and viscosity mu (mu=eta)
    fpuvOnFluidW = vec_spread3D(FatNodeW,wallRatT,hFL,fLNum,cntW); % gives a force per unit volume on the fluid of density rho and viscosity mu (mu=eta)

    fpuvOnFluid = fpuvOnFluidA + fpuvOnFluidB + fpuvOnFluidW + g0XMompV;
    [uFL] = fluid3D_IB(uFLatT,fpuvOnFluid,aDelta,dt,eta,ipIB,imIB,hFL,rho); % Step Fluid Velocity

    uFLinterpA = vec_interp3D(uFL,r3DatTA,hFL,GNA,fLNum); % interpolated force at my nodes. 
    uFLinterpB = vec_interp3D(uFL,r3DatTB,hFL,GNB,fLNum); % interpolated force at my nodes.
    uFLinterpW = vec_interp3D(uFL,wallRatT,hFL,cntW,fLNum); % interpolated force at my nodes.

    g0XuFL = g0crossFun(uFL,hFL,ipIB,imIB);
    
    omegaA = vec_interp3D(g0XuFL,r3DatTA,hFL,GNA,fLNum); % should give an x,y,z component vector of torque
    omegaB = vec_interp3D(g0XuFL,r3DatTB,hFL,GNB,fLNum); % should give an x,y,z component vector of torque
    omegaW = vec_interp3D(g0XuFL,wallRatT,hFL,cntW,fLNum); % should give an x,y,z component vector of torque

    %%
    omegaMagA = sqrt(omegaA(:,1).^2 + omegaA(:,2).^2 + omegaA(:,3).^2);
    omegaScA = omegaA./omegaMagA;

    xa(:,jj) = xatTa + uFLinterpA(:,1).*dt; 
    ya(:,jj) = yatTa + uFLinterpA(:,2).*dt;
    za(:,jj) = zatTa + uFLinterpA(:,3).*dt;    
    
    omegaMagB = sqrt(omegaB(:,1).^2 + omegaB(:,2).^2 + omegaB(:,3).^2);
    omegaScB = omegaB./omegaMagB;

    xb(:,jj) = xatTb + uFLinterpB(:,1).*dt; 
    yb(:,jj) = yatTb + uFLinterpB(:,2).*dt;
    zb(:,jj) = zatTb + uFLinterpB(:,3).*dt;

    omegaMagW = sqrt(omegaW(:,1).^2 + omegaW(:,2).^2 + omegaW(:,3).^2);
    omegaScW = omegaW./omegaMagW;

    wallR = wallRatT + uFLinterpW.*dt; 

        %% Update the d hat bases 
    [d1newA,d2newA,d3newA] = updateDHat(GNA,dt,omegaMagA,d1a,d2a,d3a,omegaScA); % cell A
    [d1newB,d2newB,d3newB] = updateDHat(GNB,dt,omegaMagB,d1b,d2b,d3b,omegaScB); % cell B
    [d1newW,d2newW,d3newW] = updateDHat(cntW,dt,omegaMagW,d1W,d2W,d3W,omegaScW); % wall   

     d1a = d1newA;d2a = d2newA;d3a = d3newA;
     d1b = d1newB;d2b = d2newB;d3b = d3newB;
     d1W = d1newW;d2W = d2newW;d3W = d3newW;
  %%  Compute half stepped d hat basis
    [d1Halfa,d2Halfa,d3Halfa] = AMatHalfRot(d1a,d2a,d3a,GNA);
    [d1Halfb,d2Halfb,d3Halfb] = AMatHalfRot(d1b,d2b,d3b,GNB);

end % end of full step rk loop 
end % end of Runge-Kutta loop
end % end of k loop

%%
delSA = r3Da(2:GNA,:) - r3Da(1:GNA-1,:);
delS2A = delSA.^2;
MagDelSA = sqrt(sum(delS2A,2));
maxMagDelSA = max(MagDelSA);
pctDelSDiffA = (maxMagDelSA./dsA-1);
% COM(jj,:) = [mean(x(:,jj)),mean(y(:,jj)),mean(z(:,jj))];
%% Side View Plotting
if GNA>1
    Om3a = dotUnit(d2a,L1A*d1a);  
else
    Om3a = 0.0;
end
if GNB>1
    Om3b = dotUnit(d2b,L1B*d1b);  
else
    Om3b = 0.0;
end


if jj == 2 || jj == 1
    Hfig = figure(1); % nothing on first iter... did pop up on 2nd
    clf(1);
end

if jj>=3
    Hfig = figure(1); % nothing on first iter... did pop up on 2nd
    clf(1);
%     Hfig.WindowState = 'minimized'; % %minimizes plot window so it doesn't interupt typing 
end
eax = 1;
eay = 1;
eaz = round(fLNum./2,0)-1;
% bact A
ra(:,1) = xa(:,jj);
ra(:,2) = ya(:,jj);
ra(:,3) = za(:,jj);
% bact B
rb(:,1) = xb(:,jj);
rb(:,2) = yb(:,jj);
rb(:,3) = zb(:,jj);

[X,Y,Z] = sphere(60);

tubeplot2(ra,Om3a,aA,32,colorvecta,colorvecta,1)
hold on
tubeplot2(rb,Om3b,aB,32,colorvectb,colorvectb,1);axis equal;
box off;axis off;

h1a = patch(surf2patch(aA(GNA).*X+ra(GNA,1),aA(GNA).*Y+ra(GNA,2),aA(GNA).*Z+ra(GNA,3),Z)); % patch last point with sphere
set(h1a,'FaceColor',colorvecta,'EdgeColor','none'); 
h2a = patch(surf2patch(aA(1).*X+ra(1,1),aA(1).*Y+ra(1,2),aA(1).*Z+ra(1,3),Z)); % patch 1st point with sphere
set(h2a,'FaceColor',colorvecta,'EdgeColor','none');box off;axis off;

h1b = patch(surf2patch(aB(GNB).*X+rb(GNB,1),aB(GNB).*Y+rb(GNB,2),aB(GNB).*Z+rb(GNB,3),Z)); % patch last point with sphere
set(h1b,'FaceColor',colorvectb,'EdgeColor','none'); 
h2b = patch(surf2patch(aB(1).*X+rb(1,1),aB(1).*Y+rb(1,2),aB(1).*Z+rb(1,3),Z)); % patch 1st point with sphere
set(h2b,'FaceColor',colorvectb,'EdgeColor','none');box off;axis off;

axis([LmxMin LmxMax LmyMin LmyMax LmzMin LmzMax]);
title(strcat('Flow in y direction, Time:',num2str(round(Time(jj),8))));
lightangle(40,-45);

% quiver fluid
uFLx = uFL(1:eax:end,1:eay:end,1:eaz:end,1); uFLy = uFL(1:eax:end,1:eay:end,1:eaz:end,2); uFLz = uFL(1:eax:end,1:eay:end,1:eaz:end,3); 
xgr = xgrid(1:eax:end,1:eay:end,1:eaz:end);ygr = ygrid(1:eax:end,1:eay:end,1:eaz:end);zgr = zgrid(1:eax:end,1:eay:end,1:eaz:end);
quiver3(xgr,ygr,zgr,uFLx,uFLy,uFLz,'k')

% Wall Plots
[my_vertices,my_faces] = makeWallPlotVars(wallR,AeffA);
patch('Vertices', my_vertices, 'Faces', my_faces, 'FaceColor', [0.6 0.6 0.6]);

[my_vertices,my_faces] = makeWallPlotVars(wallR,AeffW);
patch('Vertices', my_vertices, 'Faces', my_faces, 'FaceColor', [0.6 0.6 0.6]);

hold off
view([mySideView1 mySideView2])
set(gcf,'Color',[1 1 1])
set(Hfig, 'Position',  [50, 50, 1000, 1000]);
mov1(:,jj-1) = getframe(Hfig);
if jj == NSteps   
filename = strcat(vidtitle1,'_Fig1','.fig');
    savefig(gcf,filename)
end
%% Top View Plotting

Hfig = figure(2);  
clf(2)
tubeplot2(ra,Om3a,aA,32,colorvecta,colorvecta,1)
hold on
tubeplot2(rb,Om3b,aB,32,colorvectb,colorvectb,1);axis equal;
box off;axis off;

h1a = patch(surf2patch(aA(GNA).*X+ra(GNA,1),aA(GNA).*Y+ra(GNA,2),aA(GNA).*Z+ra(GNA,3),Z)); % patch last point with sphere
set(h1a,'FaceColor',colorvecta,'EdgeColor','none'); 
h2a = patch(surf2patch(aA(1).*X+ra(1,1),aA(1).*Y+ra(1,2),aA(1).*Z+ra(1,3),Z)); % patch 1st point with sphere
set(h2a,'FaceColor',colorvecta,'EdgeColor','none');box off;axis off;

h1b = patch(surf2patch(aB(GNB).*X+rb(GNB,1),aB(GNB).*Y+rb(GNB,2),aB(GNB).*Z+rb(GNB,3),Z)); % patch last point with sphere
set(h1b,'FaceColor',colorvectb,'EdgeColor','none'); 
h2b = patch(surf2patch(aB(1).*X+rb(1,1),aB(1).*Y+rb(1,2),aB(1).*Z+rb(1,3),Z)); % patch 1st point with sphere
set(h2b,'FaceColor',colorvectb,'EdgeColor','none');box off;axis off;

axis([xa(1,1)-2*AeffA-AeffB xb(1,1)+1.5*AeffB LmyMin LmyMax za(1,1)-7.0*max(AeffA,AeffB) zb(end,1)+7.0*max(AeffA,AeffB)]);
view([myTopView1 myTopView2])
title(strcat('IB Version, Time:',num2str(round(Time(jj),8))));
lightangle(55,55);

% plotting only the relevant fluid lines with quiver
eax = mySkipX;
eay = mySkipYZ;
eaz = mySkipYZ;
plXSt = mean(xa(:,jj) - AeffA);
findXmin = min(xgrid(xgrid(:,1,1)>plXSt));
plXEnd = mean(xb(:,jj) + 2.*AeffB);
findXmax = max(xgrid(xgrid(:,1,1)<plXEnd));
myMinX = find(xgrid(:,1,1)==findXmin);
myMaxX = find(xgrid(:,1,1)==findXmax);

for pp = 1:ceil((myMaxX-myMinX)./mySkipX)
    uFLx = uFL(myMinX+pp*eax-1,1:eay:end,1:eaz:end,1); 
    uFLy = uFL(myMinX+pp*eax-1,1:eay:end,1:eaz:end,2); 
    uFLz = uFL(myMinX+pp*eax-1,1:eay:end,1:eaz:end,3); 
    xgr = xgrid(myMinX+pp*eax-1,1:eay:end,1:eaz:end);
    ygr = ygrid(myMinX+pp*eax-1,1:eay:end,1:eaz:end);
    zgr = zgrid(myMinX+pp*eax-1,1:eay:end,1:eaz:end);
    quiver3(xgr,ygr,zgr,uFLx,uFLy,uFLz,ShowArrowHead="on",AutoScale="on",AutoScaleFactor=0.35)
end
% Wall Plot
[my_vertices,my_faces] = makeWallPlotVars(wallR,AeffW);
patch('Vertices', my_vertices, 'Faces', my_faces, 'FaceColor', [0.6 0.6 0.6]);

hold off
set(gcf,'Color',[1 1 1])
set(Hfig, 'Position',  [1050, 50, 1000, 1000]);
    
mov2(:,jj-1) = getframe(Hfig);
if jj == NSteps   
    filename = strcat(vidtitle1,'_Fig2','.fig');
    savefig(gcf,filename)
end
%% Write Movies
if mod(jj,300) == 0 || jj == NSteps
    writerObj = VideoWriter(vidtitle1);
    writerObj.FrameRate = 20; %FPS
    open(writerObj); %Open Video Writer
    for ii=1:length(mov1)
        frame = mov1(ii); %convert image to frame of movie
        writeVideo(writerObj,frame);
    end
    close(writerObj);%close the video writer


    writerObj = VideoWriter(vidtitle2);
%     writerObj.FrameRate = round((1./samplet)./10,0); %FPS 1/10th real speed
    writerObj.FrameRate = 20; %FPS
    open(writerObj); %Open Video Writer
    for ii=1:length(mov2)
        frame = mov2(ii); %convert image to frame of movie
        writeVideo(writerObj,frame);
    end
    close(writerObj);%close the video writer
end

%% Making data for output excel

maxU = max(max(max(max(uFL))));
meanU = mean(mean(mean(mean(uFL))));
meanUY = mean(mean(mean(mean(uFL(:,:,:,2)))));

displayit = ["ITER","Time","maxDelS","Effective Length (Bact A)","mean Y","max U","mean U","mean U-Y","d1by";...
    strcat(string(jj),'/',string(NSteps)),Time(jj),pctDelSDiffA,abs(za(end,jj)-za(1,jj))+2*AeffA,mean(ya(:,jj)),maxU,meanU,meanUY,d1b(1,2)];
displayit
ForceMomBOvTime(jj-1,:) = [sum(FatNodeB(:,1)),sum(FatNodeB(:,2)),sum(FatNodeB(:,3)),...
                          sum(MatNodeB(:,1)),sum(MatNodeB(:,2)),sum(MatNodeB(:,3)),...
                          sqrt(sum((sum(FatNodeB(:,:))).^2)),sqrt(sum((sum(MatNodeB(:,:))).^2))];
%% tables
if mod(jj,300) == 0 || jj == NSteps
tlist = linspace(0,TotalTime,NSteps);
tlist = tlist(2:jj);
tlist = tlist';

%% Create tables of data
Tx = table(xa); Ty = table(ya); Tz = table(za);
rLasts = zeros(jj,3);
rLasts(:,1) = xa(end,:)';
rLasts(:,2) = ya(end,:)';
rLasts(:,3) = za(end,:)';
Txyz = table(rLasts);

infodata = [string(TotalTime),maxU,meanU,eta,dt0,...
            myRadA,myRadB,phiA,phiB,yshift,zshift,...
            LA, LeffA, LB, LeffB, AeffA,...
            A1, A2, B1, B2, ...
            GNA,fLNum,myVersion,datestring];  
                 
InfoTable = table(infodata(:,1),infodata(:,2),infodata(:,3),infodata(:,4),infodata(:,5),...
                  infodata(:,6),infodata(:,7),infodata(:,8),infodata(:,9),infodata(:,10),...
                  infodata(:,11),infodata(:,12),infodata(:,13),infodata(:,14),infodata(:,15),...
                  infodata(:,16),infodata(:,17),infodata(:,18),infodata(:,19),infodata(:,20),...
                  infodata(:,21),infodata(:,22),infodata(:,23),infodata(:,24),...
                   'VariableNames', ...
                   {'TotalTime','maxU','meanU','eta','dt0',...
                    'myRadA','myRadB','phiA','phiB','yshift','zshift',...
                    'LA', 'LeffA', 'LB', 'LeffB', 'Aeff',...
                    'A1', 'A2', 'B1', 'B2',...
                    'GridNum','fLNum','myVersion','datestring'});

TForceMomB = table(MatNodeB(:,1),MatNodeB(:,2),MatNodeB(:,3),...
                   FatNodeB(:,1),FatNodeB(:,2),FatNodeB(:,3),...
                   'VariableNames', ...
                   {'MomBX','MomBY','MomBZ',...
                   'ForceBX','ForceBY','ForceBZ'});   
TForMomOvTime = table(tlist,ForceMomBOvTime(:,1),ForceMomBOvTime(:,2),ForceMomBOvTime(:,3),...
                            ForceMomBOvTime(:,4),ForceMomBOvTime(:,5),ForceMomBOvTime(:,6),...
                            ForceMomBOvTime(:,7),ForceMomBOvTime(:,8),...
                     'VariableNames', ...
                     {'Time','ForceBX','ForceBY','ForceBZ',...
                      'MomBX','MomBY','MomBZ'...
                      'TotalForce','TotMomAboutNodes'} );

%% Write tables to excel file
writetable(InfoTable,strcat(vidtitle1,'_data.xlsx'),'Sheet','Info');
writetable(TForMomOvTime,strcat(vidtitle1,'_data.xlsx'),'Sheet','F and M Over time');
writetable(TForceMomB,strcat(vidtitle1,'_data.xlsx'),'Sheet','tFinal F and M');
writetable(Tx,strcat(vidtitle1,'_data.xlsx'),'Sheet','x');
writetable(Ty,strcat(vidtitle1,'_data.xlsx'),'Sheet','y');
writetable(Tz,strcat(vidtitle1,'_data.xlsx'),'Sheet','z');
writetable(Txyz,strcat(vidtitle1,'_data.xlsx'),'Sheet','xyz');

end % end conditional write function

end % end of jj loop

end % end of function

%% Other functions
function [scalar] = dotUnit(myfunc,eMat)
        scalar = myfunc(:,1).*eMat(:,1) + myfunc(:,2).*eMat(:,2) +myfunc(:,3).*eMat(:,3);
end    
function [scalar] = dotUnitNew(myfunc,eMat)
        scalar = myfunc(:,1,:).*eMat(:,1) + myfunc(:,2,:).*eMat(:,2) + myfunc(:,3,:).*eMat(:,3);
end    
function [d1Half,d2Half,d3Half] = AMatHalfRot(d1,d2,d3,GridNum)
    ARot = zeros(3,3,GridNum-1);
    d1Half = zeros(GridNum-1,3);
    d2Half = d1Half;
    d3Half = d1Half;
    ARotHalf = ARot;    
    for kk = 1:GridNum-1
    ARot(:,:,kk) = [d1(kk+1,1).*d1(kk,1) + d2(kk+1,1).*d2(kk,1) + d3(kk+1,1).*d3(kk,1),...
                    d1(kk+1,1).*d1(kk,2) + d2(kk+1,1).*d2(kk,2) + d3(kk+1,1).*d3(kk,2),...
                    d1(kk+1,1).*d1(kk,3) + d2(kk+1,1).*d2(kk,3) + d3(kk+1,1).*d3(kk,3);...
                    d1(kk+1,2).*d1(kk,1) + d2(kk+1,2).*d2(kk,1) + d3(kk+1,2).*d3(kk,1),...
                    d1(kk+1,2).*d1(kk,2) + d2(kk+1,2).*d2(kk,2) + d3(kk+1,2).*d3(kk,2),...
                    d1(kk+1,2).*d1(kk,3) + d2(kk+1,2).*d2(kk,3) + d3(kk+1,2).*d3(kk,3);...
                    d1(kk+1,3).*d1(kk,1) + d2(kk+1,3).*d2(kk,1) + d3(kk+1,3).*d3(kk,1),...
                    d1(kk+1,3).*d1(kk,2) + d2(kk+1,3).*d2(kk,2) + d3(kk+1,3).*d3(kk,2),...
                    d1(kk+1,3).*d1(kk,3) + d2(kk+1,3).*d2(kk,3) + d3(kk+1,3).*d3(kk,3)...
                   ];
    ARotHalf(:,:,kk) = sqrtm(ARot(:,:,kk));
    d1Half(kk,:) = ARotHalf(:,:,kk)*d1(kk,:)';
    d2Half(kk,:) = ARotHalf(:,:,kk)*d2(kk,:)';
    d3Half(kk,:) = ARotHalf(:,:,kk)*d3(kk,:)';   
    end
end

function  [f1sCD,f2sCD,f3sCD] = FirstDerCD1d(f1,f2,f3,ds)
    f1sCD(:,1) = (f1(2:end)-f1(1:end-1))./ds;
    f2sCD(:,1) = (f2(2:end)-f2(1:end-1))./ds;
    f3sCD(:,1) = (f3(2:end)-f3(1:end-1))./ds;
end

function  [f1sCD,f2sCD,f3sCD] = FirstDerCD3d(f1,f2,f3,ds)
    f1sCD(:,:) = (f1(2:end,:)-f1(1:end-1,:))./ds;
    f2sCD(:,:) = (f2(2:end,:)-f2(1:end-1,:))./ds;
    f3sCD(:,:) = (f3(2:end,:)-f3(1:end-1,:))./ds;
end
function [xterm,yterm,zterm] = crossProduct(term1,term2)
    xterm = term1(:,2).*term2(:,3) - term1(:,3).*term2(:,2);
    yterm = term1(:,3).*term2(:,1) - term1(:,1).*term2(:,3);
    zterm = term1(:,1).*term2(:,2) - term1(:,2).*term2(:,1);
end
function [xterm,yterm,zterm] = crossProductNew(term13D,term2)
    xterm = term13D(:,2,:).*term2(:,3) - term13D(:,3,:).*term2(:,2);
    yterm = term13D(:,3,:).*term2(:,1) - term13D(:,1,:).*term2(:,3);
    zterm = term13D(:,1,:).*term2(:,2) - term13D(:,2,:).*term2(:,1);
end
function g0XuFL = g0crossFun(myFunc,hFL,ipIB,imIB)
            g0XuFL(:,:,:,1) = (1/(4*hFL))*(myFunc(:,ipIB,:,3) - myFunc(:,imIB,:,3) - myFunc(:,:,ipIB,2) + myFunc(:,:,imIB,2));
            g0XuFL(:,:,:,2) = (1/(4*hFL))*(myFunc(:,:,ipIB,1) - myFunc(:,:,imIB,1) - myFunc(ipIB,:,:,3) + myFunc(imIB,:,:,3));
            g0XuFL(:,:,:,3) = (1/(4*hFL))*(myFunc(ipIB,:,:,2) - myFunc(imIB,:,:,2) - myFunc(:,ipIB,:,1) + myFunc(:,imIB,:,1));
end
