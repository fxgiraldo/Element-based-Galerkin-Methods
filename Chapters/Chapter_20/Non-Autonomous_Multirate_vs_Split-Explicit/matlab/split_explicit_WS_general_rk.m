%---------------------------------------------------------------------%
%This routine advances the solution in time using Split-Explicit using 
%a simple General Order RK method.
%Written by F.X. Giraldo on 9/6/19
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qp] = split_explicit_WS_general_rk(q0,c,I,M,dt,time)

%build RK Coefficients
alphaI=zeros(I+1,1);
for i=1:I
    alphaI(i+1)=1.0/(I - i + 1);
end
cI=alphaI;

%Outer RK3 Loop
Qi=q0;
dtau=dt/M;
for i=2:I+1
    Qk=q0;
    t=time + dt*cI(i-1);
    [rhs,f_i,g_i] = rhs_function(Qi,c,t);
    
    %Inner Euler Loop
    for k=1:cI(i)*M
        t=t+dtau*(k-1);
        [rhs,f_k,g_k] = rhs_function(Qk,c,t);
        Qk=Qk + dtau*( f_k + g_i );
    end
    Qi=Qk;
end
qp=Qk;
   
   
