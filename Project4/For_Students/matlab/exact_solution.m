%---------------------------------------------------------------------%
%This function computes the Initial and Analytic Solutions.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qe,fe] = exact_solution(coord,npoin,c)

%Initialize
qe=zeros(npoin,1);
fe=zeros(npoin,1);
cc=pi*c;

%Generate Grid Points
for i=1:npoin
    x=coord(1,i);
    y=coord(2,i);
    qe(i)=sin(cc*x)*sin(cc*y);
    fe(i)=-2*cc^2*sin(cc*x)*sin(cc*y);
end %ip      


      
