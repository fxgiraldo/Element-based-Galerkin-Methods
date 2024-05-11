%----------------------------------------------------------------------%
%This subroutine constructs the flux integral contribution for the Weak Form CG/DG
%on Quadrilateral Elements for the 2D Advection Equation.
%Written by Francis X. Giraldo on 1/2001
%           Naval Postgraduate School
%           Department of Applied Mathematics
%           Monterey, CA 93943-5216
%----------------------------------------------------------------------%
function rhs = create_rhs_flux(rhs,q,u,v,face,normals,jac_face, ...
               wnq,nface,ngl,mapL,mapR,intma,flux_method,form_method)

%local arrays
ql=zeros(ngl,1);
qr=zeros(ngl,1);
ul=zeros(ngl,1);
vl=zeros(ngl,1);

flag=0;
if strcmp(flux_method,'rusanov') 
    flag=1;
end

delta=0;
if strcmp(form_method,'strong') 
    delta=1;
end

%Construct flux integral
for is=1:nface

   %Store Left Side Variables
   el=face(3,is);
   if (el ~= -6) %periodic bc
      ilocl=face(1,is);
      for l=1:ngl
          %Get Pointers
          il=mapL(1,l,ilocl);
          jl=mapL(2,l,ilocl);
          IL=intma(il,jl,el);
          
          %Left Element
          ql(l)=q(IL);
          ul(l)=u(IL);
          vl(l)=v(IL);
      end %l  

      %Store Right Side Variables
      er=face(4,is);
      if (er > 0 ) 
         ilocr=face(2,is);            
         for l=1:ngl
             %Get Pointers
             ir=mapR(1,l,ilocr);
             jr=mapR(2,l,ilocr);
             IR=intma(ir,jr,er);
             
             %Right Element
             qr(l)=q(IR);
%              ur(l)=u(IR); %Not needed since U is continuous
%              vr(l)=v(IR); %Not needed since V is continuous
         end %l
      end %if ier

      %Do Gauss-Lobatto Integration
      for l=1:ngl
          wq=wnq(l)*jac_face(l,is);
          
          %Store Normal Vectors
          nxl=normals(1,l,is);
          nyl=normals(2,l,is);
          nxr=-nxl;
          nyr=-nyl;

          %Interpolate onto Quadrature Points
          qlq_k=ql(l); qrq_k=qr(l); u_k=ul(l); v_k=vl(l);    
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
          flux_q=0.5*(flux_q - flag*diss_q);
          
          %Loop through Side Interpolation Points
          %--------------Left Side------------------%
          il=mapL(1,l,ilocl);
          jl=mapL(2,l,ilocl);
          IL=intma(il,jl,el);
          rhs(IL)=rhs(IL) - wq*(flux_q - delta*flux_ql);

          %--------------Right Side------------------%
          if (er > 0)
             ir=mapR(1,l,ilocr);
             jr=mapR(2,l,ilocr);
             IR=intma(ir,jr,er);
             rhs(IR)=rhs(IR) + wq*(flux_q + delta*flux_qr);
          end %if 
      end %l   
   end %if el      
end %is

