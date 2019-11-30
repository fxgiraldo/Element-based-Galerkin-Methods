%---------------------------------------------------------------------%
%This routine advances the solution in time using Split-Explicit using 
%a General Order RK method in Butcher tableau form following Schlegel et al. GMD 2012.
%Written by F.X. Giraldo on 9/6/19
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qp] = multirate_erk_butcher_substepping(q0,c,alphaI,betaI,I,alphaJ,betaJ,J,M,dt,time)

%Compute Coefficients
%Outer Loop
alphawI=zeros(I+1,I);
cwI=zeros(I+1,1);
cI=zeros(I,1);
for i=2:I
   cI(i)=sum(alphaI(i,:));
   alphawI(i,:)=alphaI(i,:) - alphaI(i-1,:);
   cwI(i)=cI(i) - cI(i-1);
end
for j=1:I
    alphawI(I+1,j)=betaI(j) - alphaI(I,j);
end
cwI(I+1)=1.0 - cI(I);

%Inner Loop
alphawJ=zeros(J+1,J);
cwJ=zeros(J+1,1);
cJ=zeros(J,1);
for i=2:J
   cJ(i)=sum(alphaJ(i,:));
   alphawJ(i,:)=alphaJ(i,:) - alphaJ(i-1,:);
   cwJ(i)=cJ(i) - cJ(i-1);
end
for j=1:J
    alphawJ(J+1,j)=betaJ(j) - alphaJ(J,j);
end
cwJ(J+1)=1.0 - cJ(J);

%Initialize Arrays
Q = zeros(I,1);
V = zeros(I+1,J+1);

%Outer Loop
Q(1)=q0;
for i=2:I+1
   r_i=0;
   for j=1:i-1   
       t=time+dt*cI(j);
       [rhs,f_j,g_j] = rhs_function(Q(j),c,t);
       r_i=r_i + alphawI(i,j)*g_j;
   end
   
   %Substepping Loop
    qm=Q(i-1);
    for m=1:M
        V(i,1)=qm; 
        
        %Inner Loop
        for j=2:J+1
           R_sum=0;
           for k=1:j-1
               t=time+dt*( cI(i-1) + cwI(i)*cJ(k)*m/M );
               [rhs,f_k,g_k] = rhs_function(V(i,k),c,t);
               R_sum=R_sum + alphawJ(j,k)*( r_i + cwI(i)*f_k );
           end %k
           V(i,j)=V(i,j-1) + dt/M*R_sum;
        end %j
        qm=V(i,J+1);
    end %m
    Q(i)=qm;
end %i
qp=Q(I+1);    