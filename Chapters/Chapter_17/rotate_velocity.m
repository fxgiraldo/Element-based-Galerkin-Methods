%---------------------------------------------------------------------%
%This function rotates the grid and plots it.
%Written by F.X. Giraldo on 4/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%
% INPUT LIST: ua, va, original velocity components
%             npoin are the number of global points
%             grid_rotation_angle is the grid rotation in degrees
% OUTPUT LIST:
%             ua, va, rotated velocity components
%---------------------------------------------------------------------%
function [ua,va] = rotate_velocity(ua,va,npoin,grid_rotation_angle)

%Rotate Grid
alpha=grid_rotation_angle*pi/180;

for i=1:npoin
 un=cos(alpha)*ua(i) - sin(alpha)*va(i);  
 vn=sin(alpha)*ua(i) + cos(alpha)*va(i);  
 ua(i)=un;
 va(i)=vn;
end


      


      
