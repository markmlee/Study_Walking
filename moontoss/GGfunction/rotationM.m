function RotationM = rotationM(axis,theta)
% ===================Engineered by GG============================
% Rotation matrix
% https://en.wikipedia.org/wiki/Rotation_matrix
% Input : rotation axis, angle(radian)
% output : rotation matrix
% Author: KangKyu Lee (guecom@kaist.ac.kr)
% KAIST HuboLab
% Date: 09/15/2017
% ===============================================================

% ========== header ========== %
axis = axis/norm(axis);
ux = axis(1); uy = axis(2); uz = axis(3);
% ====================================================

% ========== main ========== %
RotationM = zeros(3,3);
RotationM(1,1) = cos(theta) + ux^2*(1-cos(theta));
RotationM(1,2) = ux*uy*(1-cos(theta)) - uz*sin(theta);
RotationM(1,3) = ux*uz*(1-cos(theta)) + uy*sin(theta);
RotationM(2,1) = uy*uz*(1-cos(theta)) + uz*sin(theta);
RotationM(2,2) = cos(theta) + uy^2*(1-cos(theta));
RotationM(2,3) = uy*uz*(1-cos(theta)) + ux*sin(theta);
RotationM(3,1) = uz*ux*(1-cos(theta)) - uy*sin(theta);
RotationM(3,2) = uz*uy*(1-cos(theta)) + ux*sin(theta);
RotationM(3,3) = cos(theta) + uz^2*(1-cos(theta));
% ====================================================
end