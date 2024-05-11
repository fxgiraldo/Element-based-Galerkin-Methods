function [xgl,wgl] = cheby_gauss_eigenvalue_fxg(N)

n=N-1; %N=ngl, n=nop

%See matrix 3.243 and 3.244 in Yue's thesis
an = zeros(N,1); %n=0,1,...,nop
bn = zeros(N,1);
J = zeros(N,N);

bn(2)=0.5;
bn(3:N) = 0.25;
J(1,1)=an(1); J(1,2)=sqrt(bn(2));
for i=2:N-1
    J(i,i-1)=sqrt(bn(i));
    J(i,i)=an(i);
    J(i,i+1)=sqrt(bn(i+1));
end
J(N,N-1)=bn(N-1); J(N,N)=an(N);
xgl = eig(J);

for i=1:N
   xgr(i)=cos( (2*i-1)*pi/(2*P) ); 
   wgl(i)=pi/P;
end

l2_norm=norm(xgl-xgr,2)