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
function [qp] = split_explicit_WS_rk3_original(q0,c,M,dt,time)

%1st RK Stage
Qi=q0;
t=time;
[rhs,f_i,g_i] = rhs_function(Qi,c,t);
%M=3; %Emil
dtau=dt/M;
Qk=q0;
for k=1:M/3
    t=t+dtau*(k-1);
%     [rhs,f_k,g_k] = rhs_function(Qi,c,t); %Emil: Qk->Qi
    [rhs,f_k,g_k] = rhs_function(Qk,c,t);
    Qk=Qk + dtau*(f_k + g_i);
end
Qi=Qk;

%2nd RK Stage
%M=2; %Emil
t=time+dt/3;
[rhs,f_i,g_i] = rhs_function(Qi,c,t);
dtau=dt/M;
Qk=q0;
for k=1:M/2
    t=t+dtau*(k-1);
%    [rhs,f_k,g_k] = rhs_function(Qi,c,t); %Emil: Qk->Qi    
    [rhs,f_k,g_k] = rhs_function(Qk,c,t);
    Qk=Qk + dtau*(f_k + g_i);
end
Qi=Qk;

%3rd RK Stage
%M=1; %Emil
t=time+dt/2;
[rhs,f_i,g_i] = rhs_function(Qi,c,t);
dtau=dt/M;
Qk=q0;
for k=1:M
    t=t+dtau*(k-1);
%    [rhs,f_k,g_k] = rhs_function(Qi,c,t); %Emil: Qk->Qi
    [rhs,f_k,g_k] = rhs_function(Qk,c,t);
    Qk=Qk + dtau*(f_k + g_i);
end
qp=Qk;
   
   
