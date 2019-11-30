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
function [qp] = split_explicit_WS_rk3_rk3(q0,c,M,dt,time)

Msteps=M;
%1st RK Stage
if (M==0) 
    Msteps=3;
end
Qi=q0;
t=time;
[~,~,g_i] = rhs_function(Qi,c,t);
dtau=dt/Msteps;
Qk=q0;
Qk0=Qk;
for k=1:Msteps/3
    t0=t+dtau*(k-1);
    [~,f_k,~] = rhs_function(Qk,c,t0);
    Qk=Qk0 + dtau/3*(f_k + g_i);
    t=t0+dtau/2;
    [~,f_k,~] = rhs_function(Qk,c,t);
    Qk=Qk0 + dtau/2*(f_k + g_i);
    t=t0+dtau;
    [~,f_k,~] = rhs_function(Qk,c,t);
    Qk=Qk0 + dtau*(f_k + g_i);
    Qk0=Qk;
end
Qi=Qk; %t+dt/3 solution

%2nd RK Stage
if (M==0) 
    Msteps=2;
end
t=time+dt/3;
[~,~,g_i] = rhs_function(Qi,c,t);
dtau=dt/Msteps;
Qk=q0;
Qk0=Qk;
for k=1:Msteps/2
    t0=t+dtau*(k-1);
    [~,f_k,~] = rhs_function(Qk,c,t0);
    Qk=Qk0 + dtau/3*(f_k + g_i);
    t=t0+dtau/2;
    [~,f_k,~] = rhs_function(Qk,c,t);
    Qk=Qk0 + dtau/2*(f_k + g_i);
    t=t0+dtau;
    [~,f_k,~] = rhs_function(Qk,c,t);
    Qk=Qk0 + dtau*(f_k + g_i);
    Qk0=Qk;
end
Qi=Qk; %t+dt/2 solution

%3rd RK Stage
if (M==0) 
    Msteps=1;
end
t=time+dt/2;
[~,~,g_i] = rhs_function(Qi,c,t);
dtau=dt/Msteps;
Qk=q0;
Qk0=Qk;
for k=1:Msteps
    t0=t+dtau*(k-1);
    [~,f_k,~] = rhs_function(Qk,c,t0);
    Qk=Qk0 + dtau/3*(f_k + g_i);
    t=t0+dtau/2;
    [~,f_k,~] = rhs_function(Qk,c,t);
    Qk=Qk0 + dtau/2*(f_k + g_i);
    t=t0+dtau;
    [~,f_k,~] = rhs_function(Qk,c,t);
    Qk=Qk0 + dtau*(f_k + g_i);
    Qk0=Qk;
end
qp=Qk; %t+dt solution
   
   
