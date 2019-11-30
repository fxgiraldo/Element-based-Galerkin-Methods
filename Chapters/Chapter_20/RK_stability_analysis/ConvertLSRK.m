function [a,b,c] = ConvertLSRK(A_vec,B_vec)
% A_vec = A vector of coefficients from low storage RK
% B_vec = B vector of coefficients from low storage RK4
s = length(A_vec);
a = zeros(s,s);
b = zeros(s,1);
c = zeros(s,1);
b(s,1) = B_vec(s);
for p = s:-1:2
    b(p-1,1) = A_vec(p)*b(p,1)+B_vec(p-1);
end
for i = s:-1:1
    for j = s-1:-1:1
        if j>=i
            a(i,j) = 0;
        elseif i == j+1
            a(i,j) = B_vec(j);
        else
            a(i,j) = A_vec(j+1)*a(i,j+1)+B_vec(j);
        end
    end
end
for i = 1:s
    c(i) = sum(a(i,:),2);
end

end

