% ===============================================================
% 3d Walking Pattery Generation Kajita pg 127
% Input: Step width x,y 
% Output: Footplacement generation
% 
% Author: Moonyoung Lee (ML634@kaist.ac.kr)
% KAIST Institute Humanoid Research Lab
% Date: 09/04/2017
% ===============================================================



% ======== input ========
samples = 10;
Tperiod = 1; %time to switch CoM base

Sx = [0.0 0.3 0.3 0.3 0]; %step width x
Sy = [0.2 0.2 0.2 0.2 0.2]; %step width y


%% ======== initialze ======== %% 
numSteps = length(Sx); %change to numSteps
z = 1.0; %hubo COM height
g = 9.8;
Tc = sqrt(z/g);
C = cosh(Tperiod/Tc);
S = sinh(Tperiod/Tc);

t = 0:(Tperiod+(Tperiod)/samples)/samples:Tperiod; %time array


x0InitCond = zeros(1,numSteps); %pos IC
y0InitCond = zeros(1,numSteps);
vx0InitCond = zeros(1,numSteps); %vel IC
vy0InitCond = zeros(1,numSteps);

xPosSteps = zeros(numSteps,samples); %matrix of [row = N steps, col = x(t)]
yPosSteps = zeros(numSteps,samples); %matrix of [row = N steps, col = y(t)]
xVelSteps = zeros(numSteps,samples); %matrix of [row = N steps, col = vx(t)]
yVelSteps = zeros(numSteps,samples); %matrix of [row = N steps, col = vy(t)]

xGlobalSteps = [] ;  %to add all steps 
yGlobalSteps = []; 
vxGlobalSteps = [];
vyGlobalSteps = [];
tGlobalSteps = 0:0.1:Tperiod*numSteps-0.1;

x0InitCond(1) = 0;
y0InitCond(1) = 0;
vx0InitCond(1) = 0;
vy0InitCond(1) = 0.1;


%% p footplacement
pxFootplace = zeros(1,numSteps);
pyFootplace = zeros(1,numSteps);

for n = 2:numSteps
    pxFootplace(n) = pxFootplace(n-1) + Sx(n);
    pyFootplace(n) = pyFootplace(n-1) - (-1)^(n)* Sy(n);
end

% set origin to (0,0) for 0th step 
pyFootplace = pyFootplace + 0.2;
pyFootplace = [0 pyFootplace];
pxFootplace = [0 pxFootplace];


%% get xBar, yBar
xBar = zeros(1,numSteps-1);
yBar = zeros(1,numSteps-1);

for i = 1:numSteps
    xBar(i) = ( pxFootplace(i+1) - pxFootplace(i) ) / 2 + pxFootplace(i);
    yBar(i) = ( pyFootplace(i+1) - pyFootplace(i) ) / 2 + pyFootplace(i);
    
end







%% ======== main ========


for i=1:numSteps
    
    % ========  update init vel  ========
    vx0InitCond(i) = ( (xBar(i) - C * x0InitCond(i)) - (1-C)*pxFootplace(i) ) / (Tc*S); 
    vy0InitCond(i) = ( (yBar(i) - C * y0InitCond(i)) - (1-C)*pyFootplace(i) ) / (Tc*S); 

    
    % ======== solve x,y pos & vel ========
    
    %calculate 2d x motion
    xPosSteps(i,:) = ((x0InitCond(i) - pxFootplace(i)) * cosh(t/Tc)) + (Tc * vx0InitCond(i) * sinh(t/Tc)) + pxFootplace(i);
    %calculate 2d y motion
    yPosSteps(i,:) = ((y0InitCond(i) - pyFootplace(i)) * cosh(t/Tc)) + (Tc * vy0InitCond(i) * sinh(t/Tc)) + pyFootplace(i);
    %calculate x velocity
    xVelSteps(i,:) = ( (x0InitCond(i) - pxFootplace(i))/Tc *sinh(t/Tc) + (vx0InitCond(i)* cosh(t/Tc)) );
    %calculate y velocity
    yVelSteps(i,:) = ( (y0InitCond(i) - pyFootplace(i))/Tc *sinh(t/Tc) + (vy0InitCond(i)* cosh(t/Tc)) );
    
    
    % ======== update x0, y0Init  ========
    x0InitCond(i+1) = xPosSteps(i,samples);
    y0InitCond(i+1) = yPosSteps(i,samples);
    
    % ======== calculate terminal velocity ========
    
%     vxTerminal = xBar(i) * ( C +1 ) / (Tc*S);
%     vyTerminal = ybar* ( C -1 ) / (Tc*S);

   
end

% concat all steps 
for i = 1:numSteps
    
    %update IC of next step
    xGlobalSteps = [xGlobalSteps xPosSteps(i,:) ]; %update to add footplace
    yGlobalSteps = [yGlobalSteps yPosSteps(i,:) ];
    vxGlobalSteps = [vxGlobalSteps xVelSteps(i,:) ];
    vyGlobalSteps = [vyGlobalSteps yVelSteps(i,:) ];

end




%% plot
f1 = figure;
subplot(2,2,1)
plot(tGlobalSteps,xGlobalSteps)
xlabel('t (t)')
ylabel('x (m)')
title('X step LIPM 3D')

subplot(2,2,2)
plot(tGlobalSteps,yGlobalSteps)
xlabel('t (t)')
ylabel('y (m)')
title('Y step LIPM 3D')

subplot(2,2,3)
plot(tGlobalSteps,vxGlobalSteps)
xlabel('t (t)')
ylabel('dx/dt (m)')
title('Vel X LIPM 3D')

subplot(2,2,4)
plot(tGlobalSteps,vyGlobalSteps)
xlabel('t (t)')
ylabel('dy/dt (m)')
title('Vel Y LIPM 3D')
hold on

f2 = figure;
plot(xGlobalSteps,yGlobalSteps)
xlabel('x (m)')
ylabel('y (m)')
title('Total Steps LIPM 3D')
axis([-1 2 -0.5 0.5])
% scatter(pxFootplace,pyFootplace)
