%---------------------------------------------------------------------%
%This code computes the Legendre_Gauss_Lobatto Filter Matrix
%Written by F.X. Giraldo on 4/2000
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function f = filter_init(P,xgl,xmu)

%Constants
p=P-1;
ph=floor( (p+1)/2 );
alpha=17;
order=18;

%Initialize
leg=zeros(P,P);
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
leg_inv=inv(leg);

%Compute Weight
for i=1:P
   weight(i)=exp(  - alpha*( (i-1)/p )^order );
end
ibeg=P-4;
iend=P;
idelta=iend-ibeg;
idelta=0.05/idelta;
for i=1:ibeg-1
%    weight(i)=1;
end
for i=ibeg:iend
%    weight(i)=1.0 - (i-ibeg)*idelta;
end
weight;

%ibeg=round(2/3*P);
ibeg=1;
iend=P;
for i=1:ibeg
    weight(i)=1;
end
for i=ibeg+1:iend
    weight(i)=1.0 - ((i-ibeg)/(iend-ibeg))^2;
end
weight;

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


      
