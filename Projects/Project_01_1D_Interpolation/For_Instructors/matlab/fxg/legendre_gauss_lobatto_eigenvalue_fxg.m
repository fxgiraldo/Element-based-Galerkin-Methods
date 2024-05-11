%---------------------------------------------------------------------%
%This code computes the Legendre-Gauss-Lobatto points and weights using
%Golub-Welsch 1969 based on Yue's PhD thesis
%Written by F.X. Giraldo on 12/2023
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [xgl,wgl] = legendre_gauss_lobatto_eigenvalue_fxg(N)

n=N-1; %N=ngl, n=nop

a0=0;
an = zeros(n,1); %n=0,1,...,nop
bn = zeros(n,1); %n=1,...,nop
J = zeros(N,N);

u0=2;
a0=0;
an(:)=0;
p=1:n;
bn=p.^2 ./ (4.*p.^2-1);
bn(n)=n/(2*n-1); %called bn* in Yue's thesis and is due to fixed endpoints
J(1,1)=a0; J(1,2)=sqrt(bn(1));
for i=2:n
    J(i,i-1)=sqrt(bn(i-1));
    J(i,i)=an(i-1);
    J(i,i+1)=sqrt(bn(i));
end
J(N,N-1)=sqrt(bn(n)); J(N,N)=an(n);
[V,D] = eig(J);
for i=1:N
    xgl(i)=D(i,i);
    m1(i)=V(1,i);
end
wgl=u0*m1.^2;