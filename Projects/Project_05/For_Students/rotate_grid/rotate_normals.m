%---------------------------------------------------------------------%
%This function rotates the grid and plots it.
%Written by F.X. Giraldo on 4/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%
% INPUT LIST: coord are the coordinates
%             intma is the connectivity list
%             npoin are the number of global points
%             nelem are the number of elements
%             ngl is the number of interpolation points in an element
%             (polynomial order + 1)\
%             plot_grid is a switch to either plot or not plot
%             grid_rotation_angle is the grid rotation in degrees
% OUTPUT LIST:
%             coord are the new rotated coordinates: x=coord(:,1) and y=coord(:,2)
%---------------------------------------------------------------------%
function [nx_rotated,ny_rotated] = rotate_normals(nx,ny,nface,nq,grid_rotation_angle)

nx_rotated=zeros(nface,nq);
ny_rotated=zeros(nface,nq);

%Rotate Grid
alpha=grid_rotation_angle*pi/180;
for l=1:nface
    for i=1:nq
        x=nx(l,i);
        y=ny(l,i);
        nx_rotated(l,i)=cos(alpha)*x - sin(alpha)*y;  
        ny_rotated(l,i)=sin(alpha)*x + cos(alpha)*y;  
    end
end




      


      
