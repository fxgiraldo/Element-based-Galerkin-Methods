%---------------------------------------------------------------------%
%This routine advances the solution in time using 
%a simple General Order RK method.
%Written by F.X. Giraldo on 9/6/19
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qp] = WS_general_rk(q0,c,I,dt,time)

%build RK Coefficients
a_i=zeros(I+1,1);
for i=1:I
    a_i(i+1)=1.0/(I - i + 1);
end

qi=q0;
[rhs,f,g] = rhs_function(qi,c,time);
for i=1:I
    %Update
    qi=q0 + dt*a_i(i+1)*rhs;
    
    %Construct RHS
    [rhs,f,g] = rhs_function(qi,c,time+dt*a_i(i+1));
end %for i

%final update
qp=qi;
   
   
