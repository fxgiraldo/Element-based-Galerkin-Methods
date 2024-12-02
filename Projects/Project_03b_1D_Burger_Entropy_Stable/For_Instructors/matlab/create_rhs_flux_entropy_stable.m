%---------------------------------------------------------------------%
%This function computes the Energy-Stable and Centered fluxes in an edge-by-edge approach
%For both Strong and Weak form CGDG
%Written by F.X. Giraldo on 5/2021
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_rhs_flux_entropy_stable(rhs,qp,intma,nelem,ngl,diss,flux_method,cgdg_method)

if strcmp(cgdg_method,'strong')
    delta=1.0;
elseif strcmp(cgdg_method,'weak')
    delta=0.0;
end
         
%Integrate Flux Terms
for e=0:nelem
    L=e;
    R=e+1;
   
%     %Periodicity Flux
%     if (L < 1)
%         L=nelem;
%     end
%     if (R > nelem)
%         R=1;
%     end
    
    %Store Left State
    if (L>=1)
        IL=intma(ngl,L);
        q_l=qp(IL,1);
    else
        IL=0;
        q_l=0; %homogeneous Dirichlet BC
    end
    
    %Store Right State
    if (R<=nelem)
        IR=intma(1,R);
        q_r=qp(IR,1);
    else
        IR=0;
        q_r=0; %homogeneous Dirichlet BC
    end

    %Rusanov Flux
    %Compute Wave Speeds 
    clam=max(q_l,q_r);

    %Momentum Flux
    if strcmp(flux_method,'ES')
         %Entropy-Stable Flux
         flux=( q_l*q_l + q_l*q_r + q_r*q_r )/6.0;       
     else
         %Centered Flux
         flux=( 0.5*q_l^2 + 0.5*q_r^2 )/2.0; 
     end

    %Centered Flux
    f_l=0.5*q_l^2;
    f_r=0.5*q_r^2;
    flux_star=flux - 0.5*diss*abs(clam)*(q_r - q_l);   

    %Store RHS
    if IL > 0
        rhs(IL,1)=rhs(IL,1) - (flux_star - delta*f_l);
    end
    if IR > 0
        rhs(IR,1)=rhs(IR,1) + (flux_star - delta*f_r);
    end
end %e