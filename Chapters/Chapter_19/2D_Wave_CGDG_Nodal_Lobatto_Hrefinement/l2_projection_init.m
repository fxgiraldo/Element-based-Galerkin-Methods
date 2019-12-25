%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [P1g,P2g,P1s,P2s] = l2_projection_init(ngl,nq)

[xgl,wgl] = legendre_gauss_lobatto(ngl);
[xq,wq] = legendre_gauss_lobatto(nq);

L = lagrange_poly(xq,xgl);

%mass matrix
M = zeros(ngl);
for i=1:ngl
    for j=1:ngl
        for q=1:nq
            M(i,j) = M(i,j) + L(i,q)*L(j,q)*wq(q);
        end
    end                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
end

%integral matrices
s=0.5;
o1=-0.5;
o2=0.5;

S1g = zeros(ngl);
S2g = zeros(ngl);
S1s = zeros(ngl);
S2s = zeros(ngl);

L1 = lagrange_poly(o1+s*xq,xgl);
L2 = lagrange_poly(o2+s*xq,xgl);

for i=1:ngl
    for j=1:ngl
        for q=1:nq
           S1g(i,j) = S1g(i,j) + s*L(i,q)*L1(j,q)*wq(q);
           S2g(i,j) = S2g(i,j) + s*L(i,q)*L2(j,q)*wq(q); 
           S1s(i,j) = S1s(i,j) + L(j,q)*L1(i,q)*wq(q);
           S2s(i,j) = S2s(i,j) + L(j,q)*L2(i,q)*wq(q); 
        end
    end
end

M
S1g
S2g
S1s
S2s
inv(M)


%projection matrices

P1g = S1g*inv(M);
P2g = S2g*inv(M);
P1s = S1s*inv(M);
P2s = S2s*inv(M);

end
