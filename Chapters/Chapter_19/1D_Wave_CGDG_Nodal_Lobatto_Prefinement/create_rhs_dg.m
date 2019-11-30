%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_rhs_dg(qp,jac,nelem,ngl,wgl_matrix,psi_matrix,dpsi_matrix,u,diss,element_order,space_method)

%Initialize
rhs=zeros(ngl,nelem);
q_e=zeros(ngl,1);
wnq=zeros(ngl,1);
psi=zeros(ngl,ngl);
dpsi=zeros(ngl,ngl);

%Integrate Divergence of Flux
for e=1:nelem
   
   %Find Order of Element
   nop=element_order(e);
   n=nop+1;
   wnq(1:n)=wgl_matrix(1:n,nop);
   psi(1:n,1:n)=psi_matrix(1:n,1:n,nop);
   dpsi(1:n,1:n)=dpsi_matrix(1:n,1:n,nop);
   
   %Store Coordinates
   q_e(1:n)=qp(1:n,e);
   
   %Jacobians
   jac_e=jac(e);
   ksi_x=1/jac_e;
   
   %LGL Integration
   for l=1:n
      wq=wnq(l)*jac_e;
            
      %Form Derivative
      q_k=0;
      for j=1:n
          q_k=q_k + psi(j,l)*q_e(j)*u;
      end
      
      %Form RHS
      for i=1:n
         dhdx_i=dpsi(i,l)*ksi_x;
         rhs(i,e)=rhs(i,e) + wq*dhdx_i*q_k;
      end %i
   end %l
end %e

%No need for Fluxes if CG
if strcmp(space_method,'cg') 
    return
end

%Integrate Flux Terms
for e=1:nelem
   iel=e;
   ier=e+1;
   if (e == nelem) 
      ier=1;
   end 
   
   %Find Order of Element
   nop=element_order(e);
   n=nop+1;

   %LGL Integration
   q_l=qp(n,iel);
   q_r=qp(1,ier);
   f_l=q_l*u;
   f_r=q_r*u;
   clam=u;
   
   %Flux
   flux=0.5*( f_l + f_r - diss*abs(clam)*(q_r - q_l) );
   
   %Add to RHS
   rhs(n,iel)=rhs(ngl,iel) - flux;
   rhs(1,ier)=rhs(1,ier)     + flux;
end %ie
