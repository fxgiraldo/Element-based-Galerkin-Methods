function [P1g,P2g,P3g,P4g,P1s,P2s,P3s,P4s] = l2_projection_2D_init(ngl,nq)


[xgl,wgl] = legendre_gauss_lobatto(ngl);
[xq,wq] = legendre_gauss_lobatto(nq);

L = lagrange_poly(xq,xgl);

nngl=ngl*ngl;

%mass matrix
M = zeros(nngl);
%for n=1:nngl
for i=1:ngl
for j=1:ngl
    n=(j-1)*ngl+i;        
%     for m=1:nngl
    for k=1:ngl
    for l=1:ngl
        m=(l-1)*ngl+k;
        for p=1:nq
        for q=1:nq 
            M(n,m) = M(n,m) + L(i,p)*L(j,q)*L(k,p)*L(l,q)*wq(p)*wq(q);
        end
        end
    end
    end
end
end

M

%integral matrices
s=0.5;
o1=-0.5;
o2=0.5;

S1g = zeros(nngl);
S2g = zeros(nngl);
S3g = zeros(nngl);
S4g = zeros(nngl);
S1s = zeros(nngl);
S2s = zeros(nngl);
S3s = zeros(nngl);
S4s = zeros(nngl);

L1 = lagrange_poly(o1+s*xq,xgl);
L2 = lagrange_poly(o2+s*xq,xgl);

%for n=1:nngl
for i=1:ngl
for j=1:ngl
    n=(j-1)*ngl+i;        
%     for m=1:nngl
    for k=1:ngl
    for l=1:ngl
        m=(l-1)*ngl+k;
        for p=1:nq
        for q=1:nq 
            S1g(n,m) = S1g(n,m) + s*s*L(i,p)*L(j,q)*L1(k,p)*L1(l,q)*wq(p)*wq(q);
            S2g(n,m) = S2g(n,m) + s*s*L(i,p)*L(j,q)*L2(k,p)*L1(l,q)*wq(p)*wq(q);
            S3g(n,m) = S3g(n,m) + s*s*L(i,p)*L(j,q)*L2(k,p)*L2(l,q)*wq(p)*wq(q);
            S4g(n,m) = S4g(n,m) + s*s*L(i,p)*L(j,q)*L1(k,p)*L2(l,q)*wq(p)*wq(q);            
            S1s(n,m) = S1s(n,m) + L1(i,p)*L1(j,q)*L(k,p)*L(l,q)*wq(p)*wq(q);
            S2s(n,m) = S2s(n,m) + L2(i,p)*L1(j,q)*L(k,p)*L(l,q)*wq(p)*wq(q);
            S3s(n,m) = S3s(n,m) + L2(i,p)*L2(j,q)*L(k,p)*L(l,q)*wq(p)*wq(q);
            S4s(n,m) = S4s(n,m) + L1(i,p)*L2(j,q)*L(k,p)*L(l,q)*wq(p)*wq(q);
        end
        end
    end
    end
end
end

% for i=1:ngl
%     for j=1:ngl
%         for q=1:nq
%            S1g(i,j) = S1g(i,j) + s*L(i,q)*L1(j,q)*wq(q);
%            S2g(i,j) = S2g(i,j) + s*L(i,q)*L2(j,q)*wq(q); 
%            S1s(i,j) = S1s(i,j) + L(j,q)*L1(i,q)*wq(q);
%            S2s(i,j) = S2s(i,j) + L(j,q)*L2(i,q)*wq(q); 
%         end
%     end
% end

%projection matrices

P1g = S1g*inv(M);
P2g = S2g*inv(M);
P3g = S3g*inv(M);
P4g = S4g*inv(M);
P1s = S1s*inv(M);
P2s = S2s*inv(M);
P3s = S3s*inv(M);
P4s = S4s*inv(M);
end