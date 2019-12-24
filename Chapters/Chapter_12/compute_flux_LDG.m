%----------------------------------------------------------------------%
%This subroutine builds the FLUX term for the Weak Form LDG
%on Quadrilateral Elements.
%Written by Francis X. Giraldo on 1/2001
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%----------------------------------------------------------------------%
function rhs = compute_flux_LDG(q,psideh,nx,ny,jac_side,psi,npoin,...
               nside,ngl,nq,imapl,imapr,intma,qe,alpha,beta)

%Initialize
rhs=zeros(npoin,2);

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

%         %Interpolate Left-State onto Quadrature Points                
%         ql_k = compute_flux_LDG_qterms(ql,psi,ngl,l);         
%         %Interpolate Right-State onto Quadrature Points
%         qr_k = compute_flux_LDG_qterms(qr,psi,ngl,l);
         
        ql_k=0; qr_k=0;
        for i=1:ngl
            ql_k=ql_k + psi(i,l)*ql(i);
            qr_k=qr_k + psi(i,l)*qr(i);
        end %i 
        
        %Construct Numerical Flux
        q_mean=alpha*ql_k + beta*qr_k;
        
        %Loop through Side Interpolation Points
        for i=1:ngl

          %Left States
          h_i=psi(i,l);          
          
          %--------------Left Side------------------%
          il=imapl(ilocl,1,i);
          jl=imapl(ilocl,2,i);
          IL=intma(el,il,jl);
          rhs(IL,1)=rhs(IL,1) + wq*nxl*h_i*q_mean;
          rhs(IL,2)=rhs(IL,2) + wq*nyl*h_i*q_mean;

          %--------------Right Side------------------%
          if (er > 0)
             ir=imapr(ilocr,1,i);
             jr=imapr(ilocr,2,i);
             IR=intma(er,ir,jr);
             rhs(IR,1)=rhs(IR,1) - wq*nxl*h_i*q_mean;
             rhs(IR,2)=rhs(IR,2) - wq*nyl*h_i*q_mean;
          end %if 
        end %i
    end %l   
end %is

