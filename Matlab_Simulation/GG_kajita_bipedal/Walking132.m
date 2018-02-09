% ===================Engineered by GG============================
% Transformation matrix
% Author: KangKyu Lee (guecom@kaist.ac.kr)
% KAIST HuboLab
% Date: 10/11/2017
% ===============================================================
clear
close all
clc
addpath('GGfunction')
EngineeredbyGG
% ========== parameters ========== %
samples = 100; % display resolution
Tperiod = 1.0; % walking period
zc = 0.8;
a=1;
b=2;
Sx = [0.0 0.0 0.3 0.3 0.3 0.0]; %step width x
Sy = [0.0 0.2 0.2 0.2 0.2 0.2]; %step width y
% Sy = [0.0 0.3 0.2 0.3 0.2 0.3]; %step width y
Sth = D2R * [0 0 20 40 60 60];
% Sth = D2R * [0 0 0 0 0 0];

% ====================================================

% ========== header ========== %
numSteps = length(Sx); %change to numSteps

t = 0:Tperiod/samples:Tperiod; % ���� �𸣰ڴ�
tGlobal = 0:Tperiod/samples:(numSteps-1)*Tperiod-Tperiod/samples;

pxFootplace = zeros(1,numSteps);
pyFootplace = zeros(1,numSteps);

pxFootplacemod = zeros(1,numSteps);
pyFootplacemod = zeros(1,numSteps);

xBar = zeros(1,numSteps);
yBar = zeros(1,numSteps);


xGlobal = [];
yGlobal = [];
zGlobal = [];

vxGlobal = [];
vyGlobal = [];
vzGlobal = [];

pxstar(1) = 0;
pystar(1) = 0;
% ====================================================

%%
g=9.8;
Tc = sqrt(zc/g);
C = cosh(Tperiod/Tc);
S = sinh(Tperiod/Tc);
D = a*(C-1)^2 + b*(S/Tc)^2;

%% calculate  velocity of exchange location
x0=0;
y0=0;
x = zeros(1,samples);
y = zeros(1,samples);

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,4,[1 2 5 6])
hold on
xlim([-0.2 1.5-0.2])
ylim([-0.2 1.5-0.2])
zlim([0 1.5])
grid on
view(3)
%% ========== main ========== %
for (n=1:numSteps-1)
    % init condition�� ������ ������ �ٲ۴�
    vx0(n) = (xBar(n) - C*x0(n) - (1-C)*pxFootplacemod(n))/(Tc*S);
    vy0(n) = (yBar(n) - C*y0(n) - (1-C)*pyFootplacemod(n))/(Tc*S);
    
    for (nt=1:samples)
        %         eq4.54
        x(nt) = (x0(n) - pxFootplacemod(n)) * cosh(t(nt)/Tc) + Tc*vx0(n)*sinh(t(nt)/Tc) + pxFootplacemod(n);
        y(nt) = (y0(n) - pyFootplacemod(n)) * cosh(t(nt)/Tc) + Tc*vy0(n)*sinh(t(nt)/Tc) + pyFootplacemod(n);
        
        %         eq4.55
        vx(nt) = (x0(n) - pxFootplacemod(n))/Tc*sinh(t(nt)/Tc) + vx0(n)*cosh(t(nt)/Tc);
        vy(nt) = (y0(n) - pyFootplacemod(n))/Tc*sinh(t(nt)/Tc) + vy0(n)*cosh(t(nt)/Tc);
        
        z(nt) = zc;
        vz(nt)= 0;
    
    if(mod(nt,10)==0)
    plot3([pxFootplacemod(n) x(nt)],[pyFootplacemod(n) y(nt)],[0 zc],'k')
    end
    
    end
    x0(n+1) = x(samples);
    y0(n+1) = y(samples);
    
    % calculate footplacement
    % one step info
    pxFootplace(n+1) = pxFootplace(n) + cos(Sth(n+1)) * Sx(n+1) - sin(Sth(n+1)) * (-1)^(n+1)* Sy(n+1);
    pyFootplace(n+1) = pyFootplace(n) + sin(Sth(n+1)) * Sx(n+1) + cos(Sth(n+1)) * (-1)^(n+1)* Sy(n+1);
    % two step info
    if(n==numSteps-1)
        pxFootplace(n+2) = pxFootplace(n+1);
        pyFootplace(n+2) = pyFootplace(n+1);
    else
        pxFootplace(n+2) = pxFootplace(n+1) + cos(Sth(n+2)) * Sx(n+2) - sin(Sth(n+2)) * (-1)^(n+2)* Sy(n+2);
        pyFootplace(n+2) = pyFootplace(n+1) + sin(Sth(n+2)) * Sx(n+2) + cos(Sth(n+2)) * (-1)^(n+2)* Sy(n+2);
    end
    
    % calculate foot exchange location
    xBar(n+1) = ( pxFootplace(n+2) - pxFootplace(n+1) ) / 2 + pxFootplace(n+1);
    yBar(n+1) = ( pyFootplace(n+2) - pyFootplace(n+1) ) / 2 + pyFootplace(n+1);
    
    % eq4.59
    pxstar(n+1) = -a*(C-1)/D*(xBar(n+1) - C*(x0(n+1)));
    pystar(n+1) = -a*(C-1)/D*(yBar(n+1) - C*(y0(n+1)));
    
    % calculate foot exchange location
    pxFootplacemod(n+1) = pxFootplace(n+1) + pxstar(n+1);
    pyFootplacemod(n+1) = pyFootplace(n+1) + pystar(n+1);
           
%     vxTerminal = xBar(n) * ( C +1 ) / (Tc*S);
%     vyTerminal = yBar(n) * ( C -1 ) / (Tc*S);
    
    xGlobal = [xGlobal x];
    yGlobal = [yGlobal y];
    zGlobal = [zGlobal z];
    vxGlobal = [vxGlobal vx];
    vyGlobal = [vyGlobal vy];
    vzGlobal = [vzGlobal vz];
end
% ====================================================
plot3([pxFootplacemod(numSteps) x(end)],[pyFootplacemod(numSteps) y(end)],[0 zc],'k')


for(n=1:numSteps)
    printfoot2Dk([pxFootplace(n) pyFootplace(n) 0],Sth(n),0.2,0.15)
end
for(n=1:numSteps)
    printfoot2Dr([pxFootplacemod(n) pyFootplacemod(n) 0],Sth(n),0.2,0.15)
end
subplot(2,4,[1 2 5 6])
plot3(xGlobal,yGlobal,zGlobal)

subplot(2,4,3)
plot(tGlobal,xGlobal)
xlabel('t (t)')
ylabel('x (m)')
title('X step LIPM 3D')

subplot(2,4,4)
plot(tGlobal,yGlobal)
xlabel('t (t)')
ylabel('y (m)')
title('Y step LIPM 3D')

subplot(2,4,7)
plot(tGlobal,vxGlobal)
xlabel('t (t)')
ylabel('dx/dt (m)')
title('Vel X LIPM 3D')

subplot(2,4,8)
plot(tGlobal,vyGlobal)
xlabel('t (t)')
ylabel('dy/dt (m)')
title('Vel Y LIPM 3D')
hold on
