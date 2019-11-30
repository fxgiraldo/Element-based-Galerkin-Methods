%---------------------------------------------------------------------%
%This routine advances the solution in time using Split-Explicit using 
%a simple General Order RK method.
%Note that for M=3 in 1st stage, M=2 in 2nd, and M=1 in 3rd we need Qi in 
%RHS function instead of Qk to get 3rd order convergence.
%
%As is, we get 1st order convergence but this does not work in NUMA
%
%Written by F.X. Giraldo on 9/6/19
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qp] = split_explicit_WS_rk3_rk2_original(q0,c,M,dt,time)

%1st RK Stage
Qi=q0;
t=time;
[~,~,g_i] = rhs_function(Qi,c,t);
%M=3;
dtau=dt/M;
Qk=q0;
for k=1:M/3
    t=t+dtau*(k-1);
    [~,f_k,~] = rhs_function(Qk,c,t);
    Qk=Qk + dtau*(f_k + g_i);
end
Qi=Qk; %t+dt/3 solution

%2nd RK Stage
%M=2;
t=time+dt/3;
[~,~,g_i] = rhs_function(Qi,c,t);
dtau=dt/M;
Qk=q0;
for k=1:M/2
    t=t+dtau*(k-1);
    [~,f_k,~] = rhs_function(Qk,c,t);
    Qk=Qk + dtau*(f_k + g_i);
end
Qi=Qk; %t+dt/2 solution

%3rd RK Stage
%M=1;
t=time+dt/2;
[~,~,g_i] = rhs_function(Qi,c,t);
dtau=dt/M;
Qk=q0;
for k=1:M
    t=t+dtau*(k-1);
    [~,f_k,~] = rhs_function(Qk,c,t);
    Qk=Qk + dtau*(f_k + g_i);
end
qp=Qk; %t+dt solution
   
   
