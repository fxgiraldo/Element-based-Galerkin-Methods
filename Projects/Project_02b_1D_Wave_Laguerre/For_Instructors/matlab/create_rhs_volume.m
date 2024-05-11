%---------------------------------------------------------------------%
%This function computes the Volume Integrals for Energy-Stable flux.
%Written by F.X. Giraldo on 5/2021
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_rhs_volume(qp,u0,intma,jac,npoin,nelem,ngl,wgl,dpsi,form_method)

%Initialize
rhs=zeros(npoin,1);
weak_flag=0; strong_flag=0;
if strcmp(form_method,'weak') 
    weak_flag=1;
elseif strcmp(form_method,'strong') 
    strong_flag=1;
end

%Integrate Divergence of Flux (Strong Form)
for e=1:nelem
        
    %Jacobians
    ksi_x=1.0/jac(e);
            
    %LGL Integration
    for i=1:ngl(e)
        I=intma(i,e);
        wq=wgl(i,e)*jac(e);
        
        %Compute Divergence of Flux
        for j=1:ngl(e)
            J=intma(j,e);
            dhdx_ji=dpsi(j,i,e)*ksi_x*qp(J);
            dhdx_ij=-dpsi(i,j,e)*ksi_x*qp(I);
            dhdx=dhdx_ij*weak_flag + dhdx_ji*strong_flag;
            rhs(I)=rhs(I) - wq*u0*dhdx;
        end %j
    end %i
end %e