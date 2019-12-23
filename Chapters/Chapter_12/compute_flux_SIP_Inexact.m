%----------------------------------------------------------------------%
%This subroutine builds the FLUX vector for the Strong Form DGM-SEM
%on Quadrilateral Elements for the 2D Euler Equations.
%Written by Francis X. Giraldo on 1/2001
%           Naval Postgraduate School
%           Department of Applied Mathematics
%           Monterey, CA 93943-5216
%----------------------------------------------------------------------%
function rhs = compute_flux_SIP_Inexact(rhs,q,psideh,nx,ny,jac_side,jac,psi,dpsi,...
       nside,ngl,nq,imapl,imapr,intma,ksi_x,ksi_y,eta_x,eta_y,mu_constant,qe)

%local arrays
ql=zeros(ngl,1);
qr=zeros(ngl,1);

%Construct FVM-type Operators
for is=1:nside

    %Store Left Side Variables
    el=psideh(is,3);
    ilocl=psideh(is,1);
    for l=1:ngl
        il=imapl(ilocl,1,l);
        jl=imapl(ilocl,2,l);
        IL=intma(el,il,jl);
        ql(l)=q(IL);
    end %l  

    %Store Right Side Variables
    er=psideh(is,4);
    if (er > 0 ) 
        ilocr=psideh(is,2);            
        for l=1:ngl
             ir=imapr(ilocr,1,l);
             jr=imapr(ilocr,2,l);
             IR=intma(er,ir,jr);
             qr(l)=q(IR);
        end %l
    elseif (er == -4) %Dirichlet BC       
%         qr(:)=ql(:);
%         qr(:)=0;
        for l=1:ngl
            il=imapl(ilocl,1,l);
            jl=imapl(ilocl,2,l);
            IL=intma(el,il,jl);
            qr(l)=qe(IL);
        end %l  

    end %if ier

    %Do Gauss-Lobatto Integration
    for l=1:nq
        wq=jac_side(is,l);

        %Store Normal Vectors
        nxl=nx(is,l);
        nyl=ny(is,l);

        nxr=-nxl;
        nyr=-nyl;

        %Interpolate Left-State onto Quadrature Points
        [ql_k, ql_x, ql_y] = compute_flux_SIP_qterms_Inexact(ql,psi,dpsi,ngl,nq,...
                             ksi_x,ksi_y,eta_x,eta_y,ilocl,el,l,imapl);                  
        
        %Interpolate Right-State onto Quadrature Points
        if (er > 0)
            [qr_k, qr_x, qr_y] = compute_flux_SIP_qterms_Inexact(qr,psi,dpsi,ngl,nq,...
                                 ksi_x,ksi_y,eta_x,eta_y,ilocr,er,l,imapr);
        elseif (er == -4)
%             qr_k=0; qr_x=0; qr_y=0;
%             qr_k=ql_k; qr_x=ql_x; qr_y=ql_y;
            [qr_k, qr_x, qr_y] = compute_flux_SIP_qterms_Inexact(qr,psi,dpsi,ngl,nq,...
                                 ksi_x,ksi_y,eta_x,eta_y,ilocl,el,l,imapl);
        end
                     
        %Flux Variables
        fxl=ql_x;
        fyl=ql_y;

        fxr=qr_x;
        fyr=qr_y;

        %Normal Flux Component
        flux_ql=( nxl*fxl + nyl*fyl );
        flux_qr=( nxr*fxr + nyr*fyr );
        flux_q=flux_ql - flux_qr;

        %Dissipation Term
        mu_l=mu_constant*(ngl)*(ngl+1)/2*jac_side(is,l)/jac(el,1,l);
        mu_r=mu_constant*(ngl)*(ngl+1)/2*jac_side(is,l)/jac(el,1,l);
        mu=max(mu_l,mu_r);
        diss_q=mu*(qr_k - ql_k);

        %Construct Rusanov Flux
        flux_grad_q=0.5*(flux_q - diss_q);
        q_mean=0.5*( ql_k + qr_k );
        

        %Loop through Side Interpolation Points
        for i=1:ngl

          %Left States
          h_i=psi(i,l);          
          [psil_x, psil_y] = compute_flux_SIP_basis_Inexact(psi,dpsi,ngl,nq,...
                           ksi_x,ksi_y,eta_x,eta_y,ilocl,el,i,l,imapl);
          flux_psil=nxl*psil_x + nyl*psil_y;

          %--------------Left Side------------------%
          il=imapl(ilocl,1,i);
          jl=imapl(ilocl,2,i);
          IL=intma(el,il,jl);
          rhs(IL)=rhs(IL) + wq*h_i*flux_grad_q + wq*flux_psil*(ql_k - q_mean);

          %--------------Right Side------------------%
          if (er > 0)
             [psir_x, psir_y] = compute_flux_SIP_basis_Inexact(psi,dpsi,ngl,nq,...
                           ksi_x,ksi_y,eta_x,eta_y,ilocr,er,i,l,imapr);
             flux_psir=nxr*psir_x + nyr*psir_y;
             ir=imapr(ilocr,1,i);
             jr=imapr(ilocr,2,i);
             IR=intma(er,ir,jr);
             rhs(IR)=rhs(IR) - wq*h_i*flux_grad_q + wq*flux_psir*(qr_k - q_mean);
          end %if 
        end %i
    end %l   
end %is

