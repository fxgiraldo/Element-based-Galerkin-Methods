%---------------------------------------------------------------------%
%This function computes the Initial and Analytic Solutions.
%Written by F.X. Giraldo on 5/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qe,qe_x,qe_y,fe] = exact_solution(coord,npoin,icase)

%Initialize
qe=zeros(npoin,1);
qe_x=zeros(npoin,1);
qe_y=zeros(npoin,1);
fe=zeros(npoin,1);
c=pi;

%Generate Grid Points
for i=1:npoin
    x=coord(i,1);
    y=coord(i,2);
    if (icase == 1) %2D Solution (Homogeneous BCs)
        qe(i)=sin(c*x)*sin(c*y);
        qe_x(i)=c*cos(c*x)*sin(c*y);
        qe_y(i)=c*sin(c*x)*cos(c*y);
        fe(i)=-2*c^2*sin(c*x)*sin(c*y);
    elseif (icase == 2) %1D Solution (non-homogeneous BCs along y)
        qe(i)=sin(c*x);
        qe_x(i)=c*cos(c*x);
        qe_y(i)=0;
        fe(i)=-c^2*sin(c*x);
    elseif (icase == 3) %2D Solution (non-Homogeneous BCs)
        qe(i)=cos(c*x)*cos(c*y);
        qe_x(i)=-c*sin(c*x)*cos(c*y);
        qe_y(i)=-c*cos(c*x)*sin(c*y);
        fe(i)=-2*c^2*cos(c*x)*cos(c*y);
    end
end %ip      



      
