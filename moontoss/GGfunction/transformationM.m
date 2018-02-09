function TransformationM = transformationM(axis,theta,offset)
% ===================Engineered by GG============================
% Transformation matrix
% Input : rotation axis, angle(radian), offset
% output : rotation matrix
% Author: KangKyu Lee (guecom@kaist.ac.kr)
% KAIST HuboLab
% Date: 09/15/2017
% ===============================================================

% ========== main ========== %
TransformationM = [[rotationM(axis,theta)] [offset(1);offset(2);offset(3)];[0 0 0 1]]; 
% ====================================================

end


