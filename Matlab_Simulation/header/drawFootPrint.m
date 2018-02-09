function drawFootPrint( px, py, th)
%DRAWFOOTPRINT plots box around input foot position
%   input: foot x, foot y, rotation theta
%   output: none. Plot 
%   footprint width, height constant is set in header
%   Also use matrix rotation to draw rotated corners 

    global footWidth
    global footHeight

%   corners = [LU RU RD LD]
    LU = [px - footWidth; py + footHeight];
    RU = [px + footWidth; py + footHeight];
    RD = [px + footWidth; py - footHeight];
    LD = [px - footWidth; py - footHeight];
    
    corners = [LU RU RD LD];
    
    RotationMatrix = [cos(th) -sin(th); sin(th) cos(th)];
    
%     p' = [p] + R*[width; height]
    LU_new = [px ; py] + RotationMatrix*[-footWidth; footHeight]; 
    RU_new = [px ; py] + RotationMatrix*[footWidth; footHeight]; 
    RD_new = [px ; py] + RotationMatrix*[footWidth; -footHeight]; 
    LD_new = [px ; py] + RotationMatrix*[-footWidth; -footHeight]; 
    
%     corners_new = [LU_new RU_new RD_new LD_new];
%     X = [corners_new];
%     disp(X)
 
    plotX = [LU_new(1), RU_new(1),RD_new(1),LD_new(1)];
    plotY = [LU_new(2), RU_new(2),RD_new(2),LD_new(2)];
    
    plot3(plotX,plotY, zeros(1,length(plotX)), 'b')

end

