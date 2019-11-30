%---------------------------------------------------------------------%
%This routine advances the solution in time using Explicit RK methods 
%written in Butcher tableau form.
%Written by F.X. Giraldo on 9/6/19
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qp] = erk_butcher(q0,c,alphaI,betaI,I,dt,time)

Q = zeros(I,1);
R = zeros(I,1);

for i=1:I
   cI(i)=sum(alphaI(i,:));
end
   
Q(1)=q0;
t=time;
[rhs,f,g] = rhs_function(Q(1),c,t);
R(1)=rhs;

for i=2:I
   R_sum=0;
   for j=1:i-1
       R_sum=R_sum + alphaI(i,j)*R(j);
   end   
   Q(i)=q0 + dt*R_sum;   
   t=time+dt*cI(i);
   [rhs,f,g] = rhs_function(Q(i),c,t);
   R(i)=rhs;
end

R_sum=0;
for i=1:I
   R_sum=R_sum + betaI(i)*R(i);
end

qp=q0 + dt*R_sum;



