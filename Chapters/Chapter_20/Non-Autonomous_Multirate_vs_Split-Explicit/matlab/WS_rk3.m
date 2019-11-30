%---------------------------------------------------------------------%
%This routine advances the solution in time using Split-Explicit using 
%a simple General Order RK method.
%Written by F.X. Giraldo on 9/6/19
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qp] = WS_rk3(q0,c,dt,time)

%1st RK Stage
Q=q0;
[rhs,f,g] = rhs_function(Q,c,time);
Q=q0 + dt/3*(f + g);

%2nd RK Stage
[rhs,f,g] = rhs_function(Q,c,time + dt/3);
Q=q0 + dt/2*(f + g);

%3rd RK Stage
[rhs,f,g] = rhs_function(Q,c,time + dt/2);
Q=q0 + dt*(f + g);

%Update
qp=Q;
   
   
