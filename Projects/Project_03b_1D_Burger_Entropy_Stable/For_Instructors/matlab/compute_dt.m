%---------------------------------------------------------------------%
%This function computes the time-step to meet a given Courant_max.
%Written by F.X. Giraldo on January 19, 2024.
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [dt,u_max]=compute_dt(qe,intma,ngl,nelem,Courant_max,dx_min)

%Estimate time-step
u_max=-100000;
for e=1:nelem
   for i=1:ngl
       I=intma(i,e);
       c=qe(I,1);
       u_max=max(u_max,abs(c));
   end
end
dt=Courant_max*dx_min/u_max;