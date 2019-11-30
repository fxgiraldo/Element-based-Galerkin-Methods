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
function [coord_rotated] = rotate_grid_v2(coord,intma,npoin,nelem,ngl,plot_grid,grid_rotation_angle)

nop=ngl-1;
coord_rotated=zeros(npoin,2);

%Rotate Grid
alpha=grid_rotation_angle*pi/180;
for i=1:npoin
 xn=cos(alpha)*coord(i,1) - sin(alpha)*coord(i,2);  
 yn=sin(alpha)*coord(i,1) + cos(alpha)*coord(i,2);  
 coord_rotated(i,1)=xn;
 coord_rotated(i,2)=yn;
end


%Plot Grid
if (plot_grid == 1)
    x=zeros(5,1);
    y=zeros(5,1);
    figure;
    hold on;
    for e=1:nelem
        for j=1:ngl-1
            for i=1:ngl-1
                i1=intma(e,i,j);
                i2=intma(e,i+1,j);
                i3=intma(e,i+1,j+1);
                i4=intma(e,i,j+1);
                x(1)=coord_rotated(i1,1); y(1)=coord_rotated(i1,2);
                x(2)=coord_rotated(i2,1); y(2)=coord_rotated(i2,2);
                x(3)=coord_rotated(i3,1); y(3)=coord_rotated(i3,2);
                x(4)=coord_rotated(i4,1); y(4)=coord_rotated(i4,2);
                x(5)=coord_rotated(i1,1); y(5)=coord_rotated(i1,2);
                plot_handle=plot(x,y,'-r');
                set(plot_handle,'LineWidth',1.5);
            end
        end
        i1=intma(e,1,1);
        i2=intma(e,ngl,1);
        i3=intma(e,ngl,ngl);
        i4=intma(e,1,ngl);
        x(1)=coord(i1,1); y(1)=coord(i1,2);
        x(2)=coord(i2,1); y(2)=coord(i2,2);
        x(3)=coord(i3,1); y(3)=coord(i3,2);
        x(4)=coord(i4,1); y(4)=coord(i4,2);
        x(5)=coord(i1,1); y(5)=coord(i1,2);
        plot_handle=plot(x,y,'-b');
        set(plot_handle,'LineWidth',2);
    end
    title_text=['Rotated Grid Plot For: Ne = ' num2str(nelem) ', N = ' num2str(nop)];
    title([title_text],'FontSize',18);      

    xlabel('X','FontSize',18);
    ylabel('Y','FontSize',18);
    axis image
end 


      


      
