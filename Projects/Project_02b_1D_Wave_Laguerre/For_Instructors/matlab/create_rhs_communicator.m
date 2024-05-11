%---------------------------------------------------------------------%
%This function computes the Mass Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_rhs_communicator(rhs,Mmatrix,npoin)

for i=1:npoin
    rhs(i)=rhs(i)/Mmatrix(i);
end