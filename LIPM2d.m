% ===============================================================
% LIPM under varying initial conditions Kajita pg 110
% Input: 4 different I.C
% Output: x(t)
% 
% Author: Moonyoung Lee (ML634@kaist.ac.kr)
% KAIST Institute Humanoid Research Lab
% Date: 08/2017
% ===============================================================


% ======== input ========
x0 = -0.15;  %I.C pos
v0 = 0.46;  %I.C vel
z = 1.0; %hubo COM height
g = 9.8;
Tc = sqrt(z/g);


% ======== main ========

t = 0:0.1:1; %time array

%calculate 2d x motion
x = (x0 * cosh(t/Tc)) + (Tc * v0 * sinh(t/Tc));

%subplots
subplot(2,2,1)
xplot0 = x;
plot(t,xplot0);
title('I.C (-0.15,0.46)');

subplot(2,2,2)
x1 = -0.2;
v1 = 0.79;
xplot1 = (x1 * cosh(t/Tc)) + (Tc * v1 * sinh(t/Tc));
plot(t,xplot1);
title('I.C (-0.2,0.79)');

subplot(2,2,3)
x2 = 0.2;
v2 = -0.79;
xplot2 = (x2 * cosh(t/Tc)) + (Tc * v2 * sinh(t/Tc));
plot(t,xplot2);
title('I.C (0.2,-0.79)');

subplot(2,2,4)
x3 = 0.15;
v3 = -0.46;
xplot3 = (x3 * cosh(t/Tc)) + (Tc * v3 * sinh(t/Tc));
plot(t,xplot3);
title('I.C (0.15,-0.46)');

