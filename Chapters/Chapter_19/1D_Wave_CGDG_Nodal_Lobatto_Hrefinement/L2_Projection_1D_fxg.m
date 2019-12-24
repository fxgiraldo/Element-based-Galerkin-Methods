%----------------------------------------------------%
%This function constructs the L2 Projection
%Originally written by Michal Kopera.
%Rewritten by F.X. Giraldo based on Notes in
%Notebook #46 P. 3.31.16.6
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [P1g,P2g,P1s,P2s] = L2_Projection_1D_fxg(ngl,nq)

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

G1 = zeros(ngl);
G2 = zeros(ngl);
S1 = zeros(ngl);
S2 = zeros(ngl);

L1 = lagrange_poly(o1+s*xq,xgl);
L2 = lagrange_poly(o2+s*xq,xgl);

%Build Scatter Matrices (Refinement): From Parent to Child1 and Child2
for i=1:ngl
    for j=1:ngl
        for q=1:nq
           S1(i,j) = S1(i,j) + L(i,q)*L1(j,q)*wq(q); %Child1 
           S2(i,j) = S2(i,j) + L(i,q)*L2(j,q)*wq(q); %Child2
        end
    end
end
%Build Gather Matrices (Derefinement): From Child1 and Child2 to Parent
G1=s*S1'; %Child1
G2=s*S2'; %Child2

%projection matrices (M^{-1}S)
P1g = M\G1;
P2g = M\G2;
P1s = M\S1;
P2s = M\S2;

end
