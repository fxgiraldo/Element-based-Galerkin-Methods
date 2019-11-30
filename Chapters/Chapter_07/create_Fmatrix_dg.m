%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Maths
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function Fmatrix = create_Fmatrix_dg(inode,npoin,nelem,ngl,u,diss)

%Initialize
Fmatrix=zeros(npoin,npoin);

for e=1:nelem
    
    %Left Side
    left=e-1;
    if (e == 1)
        left=nelem;
    end
    i0=1;
    iN=ngl;
    n_left=-1;
    
    IE=inode(i0,e);
    IL=inode(iN,left);
    
    Fmatrix(IE,IE)=n_left*u/2*(1+n_left*diss);
    Fmatrix(IE,IL)=n_left*u/2*(1-n_left*diss);
    
    %Right Side
    right=e+1;
    if (e == nelem)
        right=1;
    end
    i0=1;
    iN=ngl;
    n_right=1;
    
    IE=inode(iN,e);
    IR=inode(i0,right);
    
    Fmatrix(IE,IE)=n_right*u/2*(1+n_right*diss);
    Fmatrix(IE,IR)=n_right*u/2*(1-n_right*diss);
end %e

