%---------------------------------------------------------------------%
%This code computes the Legendre_Gauss_Lobatto Filter Matrix
%Written by F.X. Giraldo on 4/2000
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [f,vdm,vdm_inv] = filter_init(P,xgl,xmu,filter_type)

%Constants
p=P-1;
ph=floor( (p+1)/2 );
alpha=17;
order=18;

%Initialize
leg=zeros(P,P);
leg2=zeros(P,P);
f=zeros(P,P);

%Compute Legendre Polynomial Matrix
for i=1:P
   x=xgl(i);
   for j=1:P
      jj=j-1;
      [L0,L0_1,L0_2]=legendre_poly(jj,x);
      leg(i,j)=L0;
   end
end

%Construct Hierarchical Modal Legendre Basis
leg2=leg;
if filter_type == 1
    for i=1:P
        x=xgl(i);
        for j=3:P
            leg2(i,j)=leg(i,j) - leg(i,j-2);
        end
    end
    leg=leg2; %Hierarchical Model Filter (no need for DSS)  
end

%Compute Inverse
leg_inv=inv(leg);

%Store Vandermonde Matrices
vdm=leg;
vdm_inv=leg_inv;

%Compute Weight
ibeg=round(2*P/3);
for i=1:ibeg
    weight(i)=1;
end
for i=ibeg+1:P
   weight(i)=exp(  - alpha*( (i-1)/p )^order );
end
weight

%Construct 1D Filter Matrix
for i=1:P
   for j=1:P
      sum=0;
      for k=1:P
         sum=sum + leg(i,k)*weight(k)*leg_inv(k,j);
      end %k
      f(i,j)=xmu*sum;
   end %j
   f(i,i)=f(i,i) + (1-xmu);
end %i	 
f

      
