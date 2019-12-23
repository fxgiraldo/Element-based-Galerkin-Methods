%----------------------------------------------------------------------%
%This subroutine builds the FLUX vector for the Strong Form DGM-SEM
%on Quadrilateral Elements for the 2D Euler Equations.
%Written by Francis X. Giraldo on 1/2001
%           Naval Postgraduate School
%           Department of Applied Mathematics
%           Monterey, CA 93943-5216
%----------------------------------------------------------------------%
function rhs = compute_flux_TensorProduct_inexact(rhs,q,u,v,face,nx,ny,jac_face,...
               nface,ngl,imapl,imapr)

%local arrays
ql=zeros(ngl,1);
qr=zeros(ngl,1);
ul=zeros(ngl,1);
vl=zeros(ngl,1);

%Construct FVM-type Operators
for is=1:nface

   %Store Left Side Variables
   iel=face(is,3);
   if (iel ~= -6) %periodic bc
      ilocl=face(is,1);
      for l=1:ngl
          %Get Pointers
          il=imapl(ilocl,1,l);
          jl=imapl(ilocl,2,l);
          %Left Element
          ql(l)=q(iel,il,jl);
          ul(l)=u(iel,il,jl);
          vl(l)=v(iel,il,jl);
      end %l  

      %Store Right Side Variables
      ier=face(is,4);
      if (ier > 0 ) 
         ilocr=face(is,2);            
         for l=1:ngl
             %Get Pointers
             ir=imapr(ilocr,1,l);
             jr=imapr(ilocr,2,l);
             %Right Element
             qr(l)=q(ier,ir,jr);
         end              %l  
      end %if ier

      %Do Gauss-Lobatto Integration
      for l=1:ngl
          wq=jac_face(is,l);
          
          %Store Normal Vectors
          nxl=nx(is,l);
          nyl=ny(is,l);

          nxr=-nxl;
          nyr=-nyl;

          %Interpolate onto Quadrature Points
          qlq_k=ql(l);
          qrq_k=qr(l);
          u_k=ul(l);
          v_k=vl(l);
             
          ul_k=u_k; vl_k=v_k;
          ur_k=u_k; vr_k=v_k;
                    
          %Compute Rusanov flux Constant
          unl=nxl*ul_k + nyl*vl_k;
          unr=nxl*ur_k + nyl*vr_k;
          claml=abs(unl);
          clamr=abs(unr);
          clam=max(claml,clamr);
          
          %Flux Variables
          fxl=qlq_k*ul_k;
          fyl=qlq_k*vl_k;

          fxr=qrq_k*ur_k;
          fyr=qrq_k*vr_k;

          %Normal Flux Component
          flux_ql=( nxl*fxl + nyl*fyl );
          flux_qr=( nxr*fxr + nyr*fyr );
          flux_q=flux_ql - flux_qr;

          %Dissipation Term
          diss_q=clam*(qrq_k - qlq_k);

          %Construct Rusanov Flux
          flux_q=0.5*(flux_q - diss_q);
            
          %--------------Left Side------------------%
          il=imapl(ilocl,1,l);
          jl=imapl(ilocl,2,l);
          rhs(iel,il,jl)=rhs(iel,il,jl) - wq*flux_q;

          %--------------Right Side------------------%
          if (ier > 0)
             ir=imapr(ilocr,1,l);
             jr=imapr(ilocr,2,l);
             rhs(ier,ir,jr)=rhs(ier,ir,jr) + wq*flux_q;
          end %if 
              
      end %l   
   end %if iel      
end %is

