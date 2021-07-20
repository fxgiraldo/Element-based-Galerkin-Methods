%---------------------------------------------------------------------%
%This function computes the Initial and Analytic Solutions.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qe,qe_x,fe] = exact_solution(coord,npoin,icase)

%Initialize
qe=zeros(npoin,1);
qe_x=zeros(npoin,1);
fe=zeros(npoin,1);
c=pi*2;

%Generate Grid Points
for ip=1:npoin
   x=coord(ip);
   if (icase == 1)
        qe(ip)=sin(c*x);
        qe_x(ip)=-c*cos(c*x);
        fe(ip)=-c^2*sin(c*x);
   end
end %ip      


      
