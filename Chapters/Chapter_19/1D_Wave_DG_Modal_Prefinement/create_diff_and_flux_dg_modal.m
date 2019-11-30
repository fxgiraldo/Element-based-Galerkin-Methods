%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_diff_and_flux_dg_modal(qp,nelem,nq,wnq,L,Ls,dL,u,diss,dx,element_order)

%Initialize
rhs=zeros(nq,nelem);
q_e=zeros(nq,1);

%Integrate Divergence of Flux
for e=1:nelem

    N=element_order(e)+1;
    %Store Coordinates
    q_e(1:N)=qp(1:N,e);

    %Jacobians
    jac=dx/2;
    ksi_x=2/dx;

    %LGL Integration
    for l=1:nq
        wq=wnq(l)*jac;

        %Form Approximation of f^e=u*q^e
        q_k=0;
        for j=1:N
            q_k=q_k + L(j,l)*q_e(j)*u;
        end

        %Form RHS
        for i=1:N
            dhdx_i=dL(i,l)*ksi_x;
            rhs(i,e)=rhs(i,e) + wq*dhdx_i*q_k;
        end %i
    end %l
end %e
%return

%Integrate Flux Terms
for e=1:nelem
    el=e;
    er=e+1;
    if (e == nelem) 
        er=1;
    end 
    Nl=element_order(el)+1;
    Nr=element_order(er)+1;
    
    %Jacobians
    jac=dx/2;
    ksi_x=2/dx;

    %Interpolate F (via q) onto element side/interface/edges
    q_l=0; q_r=0;
    
    for i=1:Nl
        q_l=q_l + Ls(i,2)*qp(i,el);
    end
    for i=1:Nr
        q_r=q_r + Ls(i,1)*qp(i,er);
    end
    f_l=q_l*u;
    f_r=q_r*u;
    clam=u;

    %Flux
    flux=0.5*( f_l + f_r - diss*abs(clam)*(q_r - q_l) );

    for i=1:Nl;
        rhs(i,el)=rhs(i,el) - Ls(i,2)*flux;
    end
    for i=1:Nr
        rhs(i,er)=rhs(i,er) + Ls(i,1)*flux;
    end
end %e
