%---------------------------------------------------------------------%
%This routine advances the solution in time using Midpoint RK2
%Written by F.X. Giraldo on 9/6/19
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qp] = rk2(q0,c,dt,time)

qi=q0;

%Evaluate at Beginning of Interval
[rhs,f,g] = rhs_function(q0,c,time);

%update to Midpoint of Interval
qi=q0 + 0.5*dt*rhs;

%update to End of Interval
[rhs,f,g] = rhs_function(qi,c,time+0.5*dt);
qi=q0 + dt*rhs;  

%update
qp=qi;
   
