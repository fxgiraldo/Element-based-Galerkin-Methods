%---------------------------------------------------------------------%
%This code computes the Laguerre-Gauss-Radau points and weights using
%Golub-Welsch 1969 based on Yue's PhD thesis
%Written by F.X. Giraldo on 12/2023
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [xgl,wgl] = Gauss_Quadrature_Eigenvalue_fxg(quad_type,N)

%quad_type=1 are Lobatto Points
%quad_type=2 are Legendre Points
%quad_type=3 are Equi-spaced Points
%quad_type=4 are Chebyshev Points
%quad_type=5 are Laguerre-Gauss-Radau Points

%Initialize
n=N-1; %N=ngl, n=nop
an = zeros(n,1); %n=0,1,...,nop
bn = zeros(n,1); %n=1,...,nop
J = zeros(N,N);
p=1:n;

%Update matrix coefficients for each Quadrature type
if (quad_type == 1) %LGL
    u0=2;
    a0=0;
    an(:)=0;
    bn=p.^2 ./ (4.*p.^2-1);
    bn(n)=n/(2*n-1); %called bn* in Yue's thesis and is due to fixed endpoints
elseif (quad_type == 2) %LG
    u0=2;
    a0=0;
    an(:)=0;
    bn=p.^2 ./ (4.*p.^2-1);
elseif (quad_type == 3) %ESP
elseif (quad_type == 4) %Chebyshev
    u0=pi;
    a0=0;
    an(:)=0;
    bn(1)=0.5;
    bn(2:n) = 0.25;
elseif (quad_type == 5) %LGR
    u0=1;
    a0=2*0+1;
    an=2.*p+1;
    an(n)=n; %called an* in Yue's thesis and is due to one fixed point (Radau rule)
    bn=p.^2;
end

%Build Jacobi matrix
J(1,1)=a0; J(1,2)=sqrt(bn(1));
for i=2:n
    J(i,i-1)=sqrt(bn(i-1));
    J(i,i)=an(i-1);
    J(i,i+1)=sqrt(bn(i));
end
J(N,N-1)=sqrt(bn(n)); J(N,N)=an(n);

%Solve Eigenvalue problem
if (quad_type <= 4)
    [V,D] = eig(J);
    for i=1:N
        xgl(i)=D(i,i);
        m1(i)=V(1,i);
    end
    wgl=u0*m1.^2; %weights fail at N >= 24 for LGR so need to do them differently
    xgl=sort(xgl);
elseif (quad_type == 5)
    xgl=sort(eig(J))';
    Lkx = laguerreL(N,xgl);
    wgl = 1./(N.*Lkx.^2);
    wgl= exp(xgl).*wgl; 
end