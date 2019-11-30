%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [Fmatrix,Fmatrix_hat,Cmatrix,Cmatrix_hat] = create_Fmatrix_hdg(intma,intma_cg,npoin,npoin_hdg,nelem,ngl,u,diss,cg_to_hdg,iperiodic_hdg)

%Initialize
Fmatrix=zeros(npoin,npoin);
Cmatrix=zeros(npoin_hdg,npoin);
Fmatrix_hat=zeros(npoin,npoin_hdg);
Cmatrix_hat=zeros(npoin_hdg,npoin_hdg);

for e=1:nelem
    
    %Left Side
    left=e-1;
    if (e == 1)
        left=nelem;
    end
    i0=1;
    iN=ngl;
    n_left=-1;
        
    %Flux Matrices
    IE=intma(i0,e);
    IL=intma(iN,left);
    Fmatrix(IE,IE)=+u*n_left*(1+n_left*diss);
    IE=intma(i0,e);
    IL=iperiodic_hdg(cg_to_hdg(intma_cg(iN,left)));
    Fmatrix_hat(IE,IL)=Fmatrix_hat(IE,IL) - u*n_left*(n_left*diss);
    
    %Conservation Condition Matrices
    IE=iperiodic_hdg(cg_to_hdg(intma_cg(i0,e)));
    IL=intma(i0,e);
    Cmatrix(IE,IL)=Cmatrix(IE,IL) + u*n_left*(1+n_left*diss);
    Cmatrix_hat(IE,IE)=Cmatrix_hat(IE,IE) - u*n_left*(n_left*diss);

    %Right Side
    right=e+1;
    if (e == nelem) %Periodic BC
        right=1;
    end
    i0=1;
    iN=ngl;
    n_right=1;
    
    %Flux Matrices
    IE=intma(iN,e);
    IR=intma(i0,right);
    Fmatrix(IE,IE)=u*n_right*(1+n_right*diss);
    IE=intma(iN,e);
    IR=iperiodic_hdg(cg_to_hdg(intma_cg(i0,right)));
    Fmatrix_hat(IE,IR)=Fmatrix_hat(IE,IR) - u*n_right*(n_right*diss);
    
    %Conservation Condition Matrices
    IE=iperiodic_hdg(cg_to_hdg(intma_cg(iN,e)));
    IR=intma(iN,e);
    Cmatrix(IE,IR)=Cmatrix(IE,IR) + u*n_right*(1+n_right*diss);
    Cmatrix_hat(IE,IE)=Cmatrix_hat(IE,IE) - u*n_right*(n_right*diss);
end %e
