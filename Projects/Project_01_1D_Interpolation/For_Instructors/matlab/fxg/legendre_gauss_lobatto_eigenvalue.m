function [xgl, wgl] = legendre_gauss_lobatto_eigenvalue(P)
% [xgl_e, wgl_e] = legendre_gauss_lobatto_eig(ngl)
% Computes the LGL points xgl_e and weights wgl_e using the 
% Eigenvalue method.  See Chapter 3 of T. Yue, 
% "Spectral element method for pricing European options and their Greeks"
%by J.F. Kelly NRL-DC

%Rearranged to match Legendre_Gauss_Lobatto

p = P - 1;                         %order
n = 1:p;
bn = n.^2 ./(4.*n.^2 -1);
bns = sqrt(bn);
bns(p) = sqrt(p./(2.*p - 1));
J = diag(bns,1) + diag(bns,-1);      %Jacobi matrix 
% size(J)
% bn
% bn(1)
% bn(p)

xgl = eig(J)';  %These are the LGL nodes
%xgl = xgl';
wgl = zeros(1,P);
for i = 1:P
    [L0,L0_1,L0_2]=legendre_poly(p,xgl(i));
    wgl(i) = 2.0./(p*(p+1).*L0^2);
end

end