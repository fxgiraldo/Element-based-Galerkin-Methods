%---------------------------------------------------------------------%
%This function computes the flux in an edge-by-edge approach (periodic).
%Written by F.X. Giraldo on 5/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_rhs_flux(rhs,qp,intma,nelem,ngl,diss,flux_method)

%Initialize
q_l=zeros(1,1);
q_r=zeros(1,1);

%Integrate Flux Terms
for e=1:nelem
    el=e;
    er=e+1;
   
    %Periodicity
    if (el < 1)
        el=nelem;
    end
    if (er > nelem)
        er=1;
    end
    
    %Store States
    il=intma(ngl,el);
    q_l(1)=qp(il,1);
    ir=intma(1,er);
    q_r(1)=qp(ir,1);

    %Rusanov Flux
    if (strcmp(flux_method,'rusanov') > 0)
        flux_star = rusanov_flux(q_l,q_r,diss);
    elseif (strcmp(flux_method,'energy') > 0)
        flux_star = energy_conserving_flux(q_l,q_r,diss);
    elseif (strcmp(flux_method,'roe') > 0)
        flux_star = roe_flux(q_l,q_r,diss);
    end

    %Store RHS
    rhs(il,1)=rhs(il,1) - flux_star(1);
    rhs(ir,1)=rhs(ir,1) + flux_star(1);
end %ie
