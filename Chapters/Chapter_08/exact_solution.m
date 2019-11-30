%---------------------------------------------------------------------%
%This function computes the Initial and Exact Solutions.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qe,qe_x,fe] = exact_solution(intma,coord,npoin,nelem,ngl,icase)

%Initialize
qe=zeros(npoin,1);
qe_x=zeros(npoin,1);
fe=zeros(npoin,1);
c=pi;

%Generate Grid Points
for e=1:nelem
    for i=1:ngl
       x=coord(i,e);
       ip=intma(i,e);
       if (icase == 1) %Homogeneous Dirichlet BCs
           c=2*pi;
           qe(ip)=sin(c*x);
           qe_x(ip)=c*cos(c*x);
           fe(ip)=-c^2*sin(c*x);
       elseif (icase == 2) %Homogeneous Dirichlet BCs
           c=pi;
           qe(ip)=cos(c*x) + 1;
           qe_x(ip)=c*sin(c*x);
           fe(ip)=-c^2*cos(c*x);
       elseif (icase == 3) %Non-Homogeneous Dirichlet BCs
           c=2*pi;
           qe(ip)=cos(c*x);
           qe_x(ip)=c*sin(c*x);
           fe(ip)=-c^2*cos(c*x);
       end
    end
end

      
