


% ======== input ========
samples = 10;
Tperiod = 1; %time to switch CoM base

Sx = [0.0 0.3 0.3 0.3 0]; %step width x
Sy = [0.2 0.2 0.2 0.2 0.2]; %step width y


% ======== initialze ========
numSteps = length(Sx); %change to numSteps
z = 1.0; %hubo COM height
g = 9.8;
Tc = sqrt(z/g);
C = cosh(Tperiod/Tc);
S = sinh(Tperiod/Tc);

% p footplacement
pxFootplace = zeros(1,numSteps);
pyFootplace = zeros(1,numSteps);

for n = 2:numSteps
    pxFootplace(n) = pxFootplace(n-1) + Sx(n);
    pyFootplace(n) = pyFootplace(n-1) - (-1)^(n)* Sy(n);
end

% set origin to (0,0) for 0th step 
pyFootplace = pyFootplace + 0.2
pyFootplace = [0 pyFootplace];
pxFootplace = [0 pxFootplace];

% walk primitive

xWalkPrimitive = zeros(1,numSteps);
yWalkPrimitive = zeros(1,numSteps);


x0InitCond = zeros(1,numSteps); %pos IC
y0InitCond = zeros(1,numSteps);
vx0InitCond = zeros(1,numSteps); %vel IC
vy0InitCond = zeros(1,numSteps);


% update primitive according to Sx, Sy
for k = 1:numSteps-1
    xWalkPrimitive(k) = Sx(k+1)/2;
    yWalkPrimitive(k) = (-1)^(k-1) * Sy(k+1)/2;
    
    % ======== calculate initial pos ========
    x0InitCond(k) = -xWalkPrimitive(k);
    y0InitCond(k) = yWalkPrimitive(k);
    
    % ======== calculate initial velocity ========
    vx0InitCond(k) = (xWalkPrimitive(k) - (x0InitCond(k)*C)) / (Tc * S);
    vy0InitCond(k) = (yWalkPrimitive(k) - (y0InitCond(k)*C)) / (Tc * S);


end



t = 0:(Tperiod+(Tperiod)/samples)/samples:Tperiod; %time array

xPosSteps = zeros(numSteps,samples); %matrix of [row = N steps, col = x(t)]
yPosSteps = zeros(numSteps,samples); %matrix of [row = N steps, col = y(t)]
xVelSteps = zeros(numSteps,samples); %matrix of [row = N steps, col = vx(t)]
yVelSteps = zeros(numSteps,samples); %matrix of [row = N steps, col = vy(t)]

xGlobalSteps = [] ; %zeros(1,samples*numSteps); %to add all steps 
yGlobalSteps = []; %zeros(1,samples*numSteps);

% ======== calculate terminal velocity ========

% 
% vxTerminal = xbar* ( C +1 ) / (Tc*S);
% vyTerminal = ybar* ( C -1 ) / (Tc*S);


% ====================== main ======================

% create walk primitive 
for i=1:numSteps-1
    
    %calculate 2d x motion
    xPosSteps(i,:) = (x0InitCond(i) * cosh(t/Tc)) + (Tc * vx0InitCond(i) * sinh(t/Tc));
    %calculate 2d y motion
    yPosSteps(i,:) = (y0InitCond(i) * cosh(t/Tc)) + (Tc * vy0InitCond(i) * sinh(t/Tc));
    %calculate x velocity
    xVelSteps(i,:) = ( x0InitCond(i)/Tc *sinh(t/Tc) + (vx0InitCond(i)* cosh(t/Tc)) );
    %calculate y velocity
    yVelSteps(i,:) = ( y0InitCond(i)/Tc *sinh(t/Tc) + (vy0InitCond(i)* cosh(t/Tc)) );

end

% concat primitives
for i = 1:numSteps-1
    
    %update IC of next step
    xGlobalSteps = [xGlobalSteps (xPosSteps(i,:) + Sx(i)*(i-1)) ]; %update to add footplace
    yGlobalSteps = [yGlobalSteps (yPosSteps(i,:)) ];
    %add to globalX,Y vs T
    
end

% ============== analy sol ==============




%============== plot ==============
f1 = figure;
subplot(3,1,1)
plot(xPosSteps(1,:),yPosSteps(1,:))
xlabel('x (m)')
ylabel('y (m)')
title('1st step LIPM 3D')

subplot(3,1,2)
plot(xPosSteps(2,:),yPosSteps(2,:))
xlabel('x (m)')
ylabel('y (m)')
title('2nd step LIPM 3D')

subplot(3,1,3)
plot(xPosSteps(3,:),yPosSteps(3,:))
xlabel('x (m)')
ylabel('y (m)')
title('3rd step LIPM 3D')


f2 = figure;
plot(xGlobalSteps,yGlobalSteps)
xlabel('x (m)')
ylabel('y (m)')
title('Total Steps LIPM 3D')
axis([-1 2 -0.5 0.5])


