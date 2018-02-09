% ===============================================================
% Header for bipedal walking 
% 
% Author: Moonyoung Lee (ML634@kaist.ac.kr)
% KAIST Institute Humanoid Research Lab
% Date: 09/06/2017
% ===============================================================


samples = 10;
Tperiod = 0.8; %time to switch CoM base

D2R = pi/180;

z = 0.8; %hubo COM height
g = 9.8;
Tc = sqrt(z/g);
C = cosh(Tperiod/Tc);
S = sinh(Tperiod/Tc);

global footWidth
global footHeight

 footWidth = 0.05; %for drawing footprint box
 footHeight = 0.025;

