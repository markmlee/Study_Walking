
% ===================Engineered by GG============================
% ============= Welcome to GG's function Library! ============== %
% Author: KangKyu Lee (guecom@kaist.ac.kr)
% KAIST HuboLab
% First Stone : 09/15/2017
% Last update : 10/12/2017
% ===============================================================
global D2R R2D 
global DPS2RPM RPM2DPS RPS2RPM RPM2RPS

global g

global Hz HzL dt dtL
%% Unit Conversion
% position conversion
D2R = pi/180;
R2D = 180/pi;
% velocity conversion
DPS2RPM = 1/360*60;
RPM2DPS = 1/DPS2RPM;
RPS2RPM = R2D*DPS2RPM;
RPM2RPS = 1/RPS2RPM;

%% Physical Parameters
g = 9.81;

%% System Parameters
% high level controller Hz
Hz = 200;
dt = 1/Hz;
% low level controller Hz
HzL = 1000;
dtL = 1/HzL;