%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_flux(rhs,qp,nelem,ngl,u,diss)

%Integrate Flux Terms
for ie=1:nelem
   iel=ie;
   ier=ie+1;
   if (ie == nelem) 
      ier=1;
   end 
   
   %LGL Integration
   q_l=qp(ngl,iel);
   q_r=qp(1,ier);
   f_l=q_l*u;
   f_r=q_r*u;
   clam=u;
   
   %Flux
   flux=0.5*( f_l + f_r - diss*abs(clam)*(q_r - q_l) );
   
   %Add to RHS
   rhs(ngl,iel)=rhs(ngl,iel) - flux;
   rhs(1,ier)=rhs(1,ier)     + flux;
end %ie
