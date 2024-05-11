%---------------------------------------------------------------------%
%This code computes the Laguerre-Gauss-Radau points and weights using
%Golub-Welsch 1969 based on Yue's PhD thesis
%Written by F.X. Giraldo on 12/2023
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [xgl,wgl] = laguerre_gauss_radau_eigenvalue_fxg(N)

n=N-1; %N=ngl, n=nop

a0=0;
an = zeros(n,1); %n=0,1,...,nop
bn = zeros(n,1); %n=1,...,nop
J = zeros(N,N);

u0=1;
a0=2*0+1;
p=1:n;
an=2.*p+1;
an(n)=n; %called an* in Yue's thesis and is due to one fixed point (Radau rule)
bn=p.^2;
J(1,1)=a0; J(1,2)=sqrt(bn(1));
for i=2:n
    J(i,i-1)=sqrt(bn(i-1));
    J(i,i)=an(i-1);
    J(i,i+1)=sqrt(bn(i));
end
J(N,N-1)=sqrt(bn(n)); J(N,N)=an(n);
xgl=sort(eig(J))';
Lkx = laguerreL(N,xgl);
wgl = 1./(N.*Lkx.^2);
wgl= exp(xgl).*wgl; 

% %Unstable Approach
% [V,D] = eig(J);
% for i=1:N
%     xgr(i)=D(i,i);
%     m1(i)=V(1,i);
% end
% wgr=u0*m1.^2; %weights fail at N >= 24
% wgr= exp(xgr).*wgr; 
% xgr=sort(xgr);
% 
% l2_norm_xgl=norm(xgl-xgr,2);
% l2_norm_wgl=norm(wgl-wgr,2);
% disp([' N = ',num2str(n),' norm(xgl-xgr) = ',num2str(l2_norm_xgl),' norm(wgl-wgr) = ',num2str(l2_norm_wgl)])
