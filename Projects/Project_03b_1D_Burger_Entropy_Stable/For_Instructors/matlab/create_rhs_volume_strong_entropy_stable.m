%---------------------------------------------------------------------%
%This function computes the Strong form volume integral for the Energy-Stable flux.
%Written by F.X. Giraldo on 5/2021
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_rhs_volume_strong_entropy_stable(qp,intma,coord,npoin,nelem,ngl,wgl,dpsi,flux_method)

%Initialize
rhs=zeros(npoin,2);
U_e=zeros(ngl,1);
x=zeros(ngl,1);
Flux=zeros(ngl,ngl);

%Integrate Divergence of Flux (Strong Form)
for e=1:nelem
    
    %Store Coordinates
    for i=1:ngl
        I=intma(i,e);
        U_e(i)=qp(I,1);
        x(i)  =coord(i,e);
    end
    
    %Jacobians
    dx=x(ngl)-x(1);
    jac=dx/2;
    ksi_x=2/dx;
    
    %Construct Flux Function
    if strcmp(flux_method,'ES')
        %Entropy-Stable Flux
        for i=1:ngl
            for j=1:ngl
                Flux(i,j)=( U_e(i)*U_e(i) + U_e(i)*U_e(j) + U_e(j)*U_e(j) )/6.0;
            end
        end
    else
        %Centered Flux
        for i=1:ngl
            for j=1:ngl
                Flux(i,j)=( 0.5*U_e(i)^2 + 0.5*U_e(j)^2 )/2.0;
            end
        end
    end
        
    %LGL Integration
    for i=1:ngl
        I=intma(i,e);
        wq=wgl(i)*jac;
        
        %Compute Divergence of Flux
        for j=1:ngl
            dhdx_ji=dpsi(j,i)*ksi_x;
            rhs(I,1)=rhs(I,1) - wq*2*dhdx_ji*Flux(i,j);
        end %j
    end %i
end %e