%----------------------------------------------------------------------%
%This subroutine constructs the flux integral contribution for the Weak Form CG/DG
%on Quadrilateral Elements for the 2D Shallow Water Equations.
%Written by Francis X. Giraldo on 11/2022
%           Naval Postgraduate School
%           Department of Applied Mathematics
%           Monterey, CA 93943-5216
%----------------------------------------------------------------------%
function rhs = create_rhs_flux(rhs,q,face,normals,jac_face, ...
               wnq,nface,ngl,mapL,mapR,intma,space_method,flux_method,form_method)

%local arrays
rl=zeros(ngl,1);
ul=zeros(ngl,1);
vl=zeros(ngl,1);
rr=zeros(ngl,1);
ur=zeros(ngl,1);
vr=zeros(ngl,1);

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
   er=face(4,is);
   if (strcmp(space_method,'dg') || (strcmp(space_method,'cg') && er < 0) )
   if (el ~= -6) %periodic bc
      ilocl=face(1,is);
      for l=1:ngl
          %Get Pointers
          il=mapL(1,l,ilocl);
          jl=mapL(2,l,ilocl);
          IL=intma(il,jl,el);
          
          %Left Element
          rl(l)=q(1,IL);
          ul(l)=q(2,IL)/rl(l);
          vl(l)=q(3,IL)/rl(l);
      end %l  

      %Store Right Side Variables
      if (er > 0 ) 
         ilocr=face(2,is);            
         for l=1:ngl
             %Get Pointers
             ir=mapR(1,l,ilocr);
             jr=mapR(2,l,ilocr);
             IR=intma(ir,jr,er);
             
             %Right Element
             rr(l)=q(1,IR);
             ur(l)=q(2,IR)/rr(l); 
             vr(l)=q(3,IR)/rr(l);
         end %l
      elseif (er <= 0 ) 
         rr=rl; ur=ul; vr=vl;
         ilocr=face(2,is);            
         for l=1:ngl
             nxl=normals(1,l,is);
             nyl=normals(2,l,is);
             unl=nxl*ul(l) + nyl*vl(l);
             rr(l)=rl(l);
             ur(l)=ur(l) - 2*unl*nxl;
             vr(l)=vr(l) - 2*unl*nyl;
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
          rl_k=rl(l);  ul_k=ul(l); vl_k=vl(l);
          rr_k=rr(l);  ur_k=ur(l); vr_k=vr(l);    
                    
          %Compute Rusanov flux Constant
          unl=nxl*ul_k + nyl*vl_k;
          unr=nxl*ur_k + nyl*vr_k;
          claml=abs(unl);
          clamr=abs(unr);
          clam=max(claml,clamr);
          
          %Flux Variables
          flux_l=zeros(3,2);
          flux_l(1,1)=rl_k*ul_k; 
          flux_l(1,2)=rl_k*vl_k;
          flux_l(2,1)=rl_k*ul_k*ul_k + 0.5*rl_k^2; 
          flux_l(2,2)=rl_k*ul_k*vl_k;
          flux_l(3,1)=rl_k*vl_k*ul_k; 
          flux_l(3,2)=rl_k*vl_k*vl_k + 0.5*rl_k^2;

          flux_r=zeros(3,2);
          flux_r(1,1)=rr_k*ur_k; 
          flux_r(1,2)=rr_k*vr_k;
          flux_r(2,1)=rr_k*ur_k*ur_k + 0.5*rr_k^2; 
          flux_r(2,2)=rr_k*ur_k*vr_k;
          flux_r(3,1)=rr_k*vr_k*ur_k; 
          flux_r(3,2)=rr_k*vr_k*vr_k + 0.5*rr_k^2;

          %Normal Flux Component
          flux_left=zeros(3,1); flux_right=zeros(3,1); flux=zeros(3,1);
          for m=1:3
            flux_left(m) =( nxl*flux_l(m,1) + nyl*flux_l(m,2) );
            flux_right(m)=( nxr*flux_r(m,1) + nyr*flux_r(m,2) );
            flux(m)=flux_left(m) - flux_right(m);
          end

          %Dissipation Term
          diss=zeros(3,1);
          diss(1)=clam*(rr_k - rl_k);
          diss(2)=clam*(rr_k*ur_k - rl_k*ul_k);
          diss(3)=clam*(rr_k*vr_k - rl_k*vl_k);

          %Construct Rusanov Flux
          flux_star=zeros(3,1);
          for m=1:3
            flux_star(m)=0.5*(flux(m) - flag*diss(m));
          end
          
          %Loop through Side Interpolation Points
          %--------------Left Side------------------%
          il=mapL(1,l,ilocl);
          jl=mapL(2,l,ilocl);
          IL=intma(il,jl,el);
          for m=1:3
            rhs(m,IL)=rhs(m,IL) - wq*(flux_star(m) - delta*flux_left(m));
          end

          %--------------Right Side------------------%
          if (er > 0)
             ir=mapR(1,l,ilocr);
             jr=mapR(2,l,ilocr);
             IR=intma(ir,jr,er);
             for m=1:3
                rhs(m,IR)=rhs(m,IR) + wq*(flux_star(m) + delta*flux_right(m));
             end
          end %if 
      end %l   
   end %if el     
   end %if space_method
end %is

