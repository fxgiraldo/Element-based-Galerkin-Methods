%---------------------------------------------------------------------%
%This function computes the Initial and Exact Solutions.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qe,u0] = exact_solution(intma,coord,npoin,nelem,nelem_LGL,ngl)

%Initialize
qe=zeros(npoin,1);

%Set some constants
xmin=coord(1,1);
xmax=coord(ngl(nelem_LGL),nelem_LGL);
%xm=0.5*(xmax + xmin);
xm=0.5*(xmax + xmin);
xc=xm;
u0=2;

%Generate Solution
for e=1:nelem
    for i=1:ngl(e)
        x=coord(i,e);
        ip=intma(i,e);
        %qe(ip)=exp(-0.1*(x-xc)^2);
        %qe(ip)=exp(-100*(x-xc)^2);
        qe(ip)=exp(-16*(x-xc)^2);
    end %i
end %e  


      
