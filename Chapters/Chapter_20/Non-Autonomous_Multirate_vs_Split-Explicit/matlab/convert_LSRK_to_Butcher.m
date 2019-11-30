function [a,b,c] = convert_LSRK_to_Butcher(A_vec,B_vec)

%Initialize Butcher Arrays
s=length(A_vec);
a = zeros(s,s);
b = zeros(s,1);
c = zeros(s,1);

%Build the b-vector
b(s,1) = B_vec(s);
for i = s:-1:2
    b(i-1,1) = A_vec(i)*b(i,1)+B_vec(i-1);
end

%Build the A-matrix
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

%Build the c-vector
for i = 1:s
    c(i) =sum(a(i,:));
end

end
    

