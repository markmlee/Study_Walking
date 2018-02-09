function Printfoot2D = printfoot2Dk(offset,theta,fx,fy)
% ===================Engineered by GG============================
% plot foot in 2D
% Input : foot position, foot angle(radian), foot dimension x, foot dimension y
% output : plot foot in 2D
% Author: KangKyu Lee (guecom@kaist.ac.kr)
% KAIST HuboLab
% Date: 09/15/2017
% ===============================================================

% ========== header ========== %
LU = [-fx/2;fy/2;0;1];
RU = [fx/2;fy/2;0;1];
RD = [fx/2;-fy/2;0;1];
LD = [-fx/2;-fy/2;0;1];
% ====================================================

% ========== main ========== %
RotationM = rotationM([0 0 1]',theta);
TransformationM = transformationM([0 0 1]',theta,offset);
LUR = TransformationM*LU;
RUR = TransformationM*RU;
RDR = TransformationM*RD;
LDR = TransformationM*LD;
footset = [LUR RUR RDR LDR LUR];

x = footset(1,:);
y = footset(2,:);

plot(x,y,'r')
% ====================================================

end