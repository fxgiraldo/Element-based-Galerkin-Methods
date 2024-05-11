function f = filter_init(xgr,xmu)
% f = filter_init(P,xgr,xmu)
% Creates a vandeven low pass filter for an infinite element
% James F. Kelly
% 25 Feb 2023

%Constants
P = length(xgr);
p=P-1;
%cutoff =  2/3;
cutoff = 1/2;
erf_order = 12;   
%Initialize
lag=zeros(P,P);
f=zeros(P,P);

%Compute a scaled Laguerre Polynomial Matrix
for i=1:P
   x=xgr(i);
   for j=1:P
      jj=j-1;
      %[L0,L0_1,L0_2]=legendre_poly(jj,x);
      lag(i,j)= scaled_laguerre1(jj,x);
   end
end
lag_inv=inv(lag);

%Compute Weight
weight = vandeven_modal(erf_order,p,cutoff,p);   %Compte the Boyd-Vandeven Transfer Function 

%Construct 1D Filter Matrix
for i=1:P
   for j=1:P
      sum=0;
      for k=1:P
         sum=sum + lag(i,k)*weight(k)*lag_inv(k,j);
      end %k
      f(i,j)=xmu*sum;
   end %j
   f(i,i)=f(i,i) + (1-xmu);
end %i	 

end

      
