%---------------------------------------------------------------------%
%This code computes the Chebyshev-Gauss points and weights using
%Golub-Welsch 1969 based on Yue's PhD thesis
%which are the roots of the Lobatto Polynomials.
%Written by F.X. Giraldo on 12/2023
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [xgl,wgl] = chebyshev_gauss_eigenvalue_fxg(N)

n=N-1; %N=ngl, n=nop

an = zeros(n,1); %n=0,1,...,nop -> an(1)=an(0)
bn = zeros(n,1); %n=1,...,nop -> bn(1)=bn(1)
J = zeros(N,N);

u0=pi;
a0=0;
an(:)=0;
bn(1)=0.5;
bn(2:n) = 0.25;
J(1,1)=an(1); J(1,2)=sqrt(bn(1));
for i=2:N-1
    J(i,i-1)=sqrt(bn(i-1));
    J(i,i)=an(i-1);
    J(i,i+1)=sqrt(bn(i));
end
J(N,N-1)=sqrt(bn(n)); J(N,N)=an(N);
[V,D] = eig(J);
for i=1:N
    xgl(i)=D(i,i);
    m1(i)=V(1,i);
end
wgl=pi*m1.^2;

% %Verify the Results
% for i=1:N
%    xgr(i)=cos( (2*i-1)*pi/(2*N) ); 
%    wgr(i)=pi/N;
% end
% xgr=sort(xgr);
% l2_norm_xgl=norm(xgl-xgr,2);
% l2_norm_wgl=norm(wgl-wgr,2);
% disp([' N = ',num2str(n),' norm(xgl-xgr) = ',num2str(l2_norm_xgl),' norm(wgl-wgr) = ',num2str(l2_norm_wgl)])
