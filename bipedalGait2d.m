% ===============================================================
% 2D walking pattern generation w/ LIPM Kajita pg 107
% Eqtn used (4.5)
% Input: I.C (pos, vel)
% Output: 2D animation of x(t) of LIPM for 2 steps
% 
% Author: Moonyoung Lee (ML634@kaist.ac.kr)
% KAIST Institute Humanoid Research Lab
% Date: 08/2017
% ===============================================================


x0InitCond = zeros(5,1);
v0InitCond = zeros(5,1);


z = 1.0; %hubo COM height
g = 9.8;
Tc = sqrt(z/g);

% ======== input ========
x0InitCond(1) = -0.2;  %I.C pos
v0InitCond(1) = 0.79;  %I.C vel
stepSize = 0.3; %how far to move LIMP base 
Tperiod = 0.5; %time to switch CoM base

% ======== main ========
N = 10;
t = 0:Tperiod/N:Tperiod; %time array

%calculate 2d x motion
xFirstStep = (x0InitCond(1) * cosh(t/Tc)) + (Tc * v0InitCond(1) * sinh(t/Tc));
vFirstStep = (x0InitCond(1)/Tc * sinh(t/Tc)) + (v0InitCond(1) * cosh(t/Tc));

%update IC
x0InitCond(2) = xFirstStep(N);  %I.C pos 2ndStep
v0InitCond(2) = vFirstStep(N);  %I.C vel 2ndStep

%add to global X vs T
xGlobal = xFirstStep;

% === next step ===

%calculate 2d x motion 
xSecondStep = ( (stepSize - x0InitCond(2)) * cosh(t/Tc)) + (Tc * v0InitCond(2) * sinh(t/Tc));
vSecondStep = (x0InitCond(2)/Tc * sinh(t/Tc)) + (v0InitCond(2) * cosh(t/Tc));


%add to global X vs T
stepOffset = xFirstStep(N) - xSecondStep(1);
xGlobalSecond = xSecondStep + stepOffset;
xGlobal = [xGlobal xGlobalSecond];


%plot
f1 = figure;


subplot(3,1,1)
plot(t,xFirstStep);
xlabel('time (t)')
ylabel('x distance (m)')
title('1st step LIPM')

subplot(3,1,2)
plot(t,xSecondStep);
xlabel('time (t)')
ylabel('x distance (m)')
title('2nd step LIPM')

subplot(3,1,3)
xTotal = [xFirstStep  xSecondStep];
plot(xTotal);
xlabel('time (t)')
ylabel('x distance (m)')
title('Combined LIPM')

f2 = figure;
hold on
xlabel('x distance (m)')
ylabel('y height (m)')
title('LIPM movement')
%animation plot
for i=1:length(xGlobal)
   pause (0.1);
   plot(xGlobal(i),z,'or','MarkerSize', 5, 'MarkerFaceColor', 'r')
   axis([-1 3 0 2])

%    delete(linePlot);
   
   if i <= (length(xGlobal))/2
       linePlot = plot([0,xGlobal(i)],[0,z],'b');
    
   else
       linePlot = plot([stepSize,xGlobal(i)],[0,z],'b');
   end 
   
   
   
end

