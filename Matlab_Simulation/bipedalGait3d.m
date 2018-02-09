


% ======== input ========
samples = 10;
Tperiod = 1; %time to switch CoM base

Sx = [0.3 0.3 0.3]; %step width x
Sy = [0.2 0.2 0.2]; %step width y


% ======== initialze ========
numSteps = length(Sx); %change to numSteps
z = 1.0; %hubo COM height
g = 9.8;
Tc = sqrt(z/g);
C = cosh(Tperiod/Tc);
S = sinh(Tperiod/Tc);

% p footplacement


% walk primitive

xWalkPrimitive = zeros(1,numSteps);
yWalkPrimitive = zeros(1,numSteps);


x0InitCond = zeros(numSteps,1); %pos IC
y0InitCond = zeros(numSteps,1);
vx0InitCond = zeros(numSteps,1); %vel IC
vy0InitCond = zeros(numSteps,1);


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



t = 0:Tperiod/samples:Tperiod-(Tperiod/samples); %time array

xPosSteps = zeros(numSteps,samples); %matrix of [row = N steps, col = x(t)]
yPosSteps = zeros(numSteps,samples); %matrix of [row = N steps, col = y(t)]
xVelSteps = zeros(numSteps,samples); %matrix of [row = N steps, col = vx(t)]
yVelSteps = zeros(numSteps,samples); %matrix of [row = N steps, col = vy(t)]



% ======== calculate terminal velocity ========

% 
% vxTerminal = xbar* ( C +1 ) / (Tc*S);
% vyTerminal = ybar* ( C -1 ) / (Tc*S);


% ======== main ========

for i=1:numSteps
    
    %calculate 2d x motion
    xPosSteps(i,:) = (x0InitCond(i) * cosh(t/Tc)) + (Tc * vx0InitCond(i) * sinh(t/Tc));
    %calculate 2d y motion
    yPosSteps(i,:) = (y0InitCond(i) * cosh(t/Tc)) + (Tc * vy0InitCond(i) * sinh(t/Tc));
    %calculate x velocity
    xVelSteps(i,:) = ( x0InitCond(i)/Tc *sinh(t/Tc) + (vx0InitCond(i)* cosh(t/Tc)) );
    %calculate y velocity
    yVelSteps(i,:) = ( y0InitCond(i)/Tc *sinh(t/Tc) + (vy0InitCond(i)* cosh(t/Tc)) );
    
    %update IC of next step
    
    %add to globalX,Y vs T
    
end


%plot
f1 = figure;
plot(xPosSteps(1,:),yPosSteps(1,:))
xlabel('x (m)')
ylabel('y (m)')
title('1st step LIPM 3D')

f2 = figure;
subplot(2,1,1)
plot(t, xVelSteps(1,:))
xlabel('time (t)')
ylabel('dx/dt (m/s)')
title('1st step LIPM 3D')

subplot(2,1,2)
plot(t, yVelSteps(1,:))
xlabel('time (t)')
ylabel('dy/dt (m/s)')
title('1st step LIPM 3D')



