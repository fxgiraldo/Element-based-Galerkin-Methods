%---------------------------------------------------------------------%
%This function computes the Initial and Analytic Solutions.
%Written by F.X. Giraldo on 5/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qe,qe_x,qe_y,qe_xx,qe_yy,fe] = exact_solution_test_derivatives(coord,npoin,icase)

%Initialize
qe=zeros(npoin,1);
qe_x=zeros(npoin,1);
qe_y=zeros(npoin,1);
qe_xx=zeros(npoin,1);
qe_yy=zeros(npoin,1);
fe=zeros(npoin,1);
c=pi;

%Generate Grid Points
for i=1:npoin
   x=coord(i,1);
   y=coord(i,2);

%    qe(i)  =sin(c*x);
%    qe_x(i)=c*cos(c*x);
%    qe_xx(i)=-c^2*sin(c*x);
%    fe(i)=-c^2*sin(c*x);

%    qe(i)  =sin(c*y);
%    qe_y(i)=c*cos(c*y);
%    qe_yy(i)=-c^2*sin(c*y);
%    fe(i)=-c^2*sin(c*y);

   
%            qe(i)=c*cos(c*x);
%            qe_x(i)=-c^2*sin(c*x);
%            qe_xx(i)=-c^3*cos(c*x);
%            fe(i)=qe_xx(i);

   qe(i)  =sin(c*x)*sin(c*y);
   qe_x(i)=c*cos(c*x)*sin(c*y);
   qe_y(i)=c*sin(c*x)*cos(c*y);
   qe_xx(i)=-c^2*sin(c*x)*sin(c*y);
   qe_yy(i)=-c^2*sin(c*x)*sin(c*y);
   fe(i)=-2*c^2*sin(c*x)*sin(c*y);
end %ip      



      
