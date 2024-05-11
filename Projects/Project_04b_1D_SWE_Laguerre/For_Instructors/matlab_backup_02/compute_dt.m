%---------------------------------------------------------------------%
%This function computes the time-step to meet a given Courant_max.
%Written by F.X. Giraldo on January 19, 2024.
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function dt=compute_dt(qe,qb,intma,ngl,nelem,icase,Courant_max,dx_min,gravity)

%Estimate time-step
u_max=-100000;
for e=1:nelem
   for i=1:ngl(e)
       I=intma(i,e);
       h=qe(I,1)+qb(I);
       u=qe(I,2)/h;
       c=u + sqrt(gravity*h);
       if (icase == -1)
           c=qe(I,2);
       end
       u_max=max(u_max,c);
   end
end
dt=Courant_max*dx_min/u_max;