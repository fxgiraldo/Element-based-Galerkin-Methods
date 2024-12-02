%---------------------------------------------------------------------%
%This function computes the Initial and Exact Solutions.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qe] = exact_solution_dg(intma,coord,npoin,nelem,ngl,time,icase)

%Initialize
qe=zeros(npoin,2);
qb=zeros(npoin,1);

%Set some constants
xmin=min(min(coord));
xmax=max(max(coord));
xm=0.5*(xmax + xmin);
xc=xm;
xl=xmax-xmin;

%Generate Grid Points
for e=1:nelem
    for i=1:ngl
        x=coord(i,e);
        ip=intma(i,e);
        if (icase == 1) %Gaussian IC with flat bottom
            qe(ip,1)=sin(pi*(x-time)) + 0.01;
            qe(ip,2)=0;
        end
    end %i
end %e  


      
