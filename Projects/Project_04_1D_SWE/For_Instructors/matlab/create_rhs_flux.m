%---------------------------------------------------------------------%
%This function computes the RHS Flux contribution for the 1D SWE.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_rhs_flux(rhs,qp,qb,intma,nelem,ngl,diss,gravity,delta_nl)

%Initialize
q_l=zeros(3,1);
q_r=zeros(3,1);

%Integrate Flux Terms
for e=0:nelem
   el=e;
   er=e+1;
   
   %LGL Integration
   if (el >= 1)
       il=intma(ngl,el);
       hs_l=qp(il,1);
       U_l =qp(il,2);
       hb_l=qb(il);
       h_l =hs_l + hb_l;
       u_l =U_l/(h_l);
   end
     
   if (er <= nelem)   
       ir=intma(1,er);
       hs_r=qp(ir,1);
       U_r =qp(ir,2);
       hb_r=qb(ir);
       h_r =hs_r + hb_r;
       u_r =U_r/(h_r);
   end
   
   %Left Boundary Conditions
   if (el == 0)
        h_l =h_r;
        hs_l=hs_r;
        U_l =-U_r; %hard wall boundaries
        hb_l=hb_r;
        u_l =U_l/(h_l);   
   end
 
   %Right Boundary Conditions
   if (er == nelem+1)
        h_r =h_l;
        hs_r=hs_l;
        U_r =-U_l; %hard wall boundaries
        hb_r=hb_l;
        u_r =U_r/(h_r);   
   end

   %Store States
   q_l(1)=hs_l;
   q_l(2)=hb_l;
   q_l(3)=u_l;
   q_r(1)=hs_r;
   q_r(2)=hb_r;
   q_r(3)=u_r;

   %Rusanov Flux
   flux = rusanov_flux(q_l,q_r,diss,gravity,delta_nl);
   
   %Add to RHS to Left
   if (el >= 1)
    rhs(il,1)=rhs(il,1) - flux(1);
    rhs(il,2)=rhs(il,2) - flux(2);
   end
   
   %Add to RHS to Right
   if (er <= nelem)
    rhs(ir,1)=rhs(ir,1) + flux(1);
    rhs(ir,2)=rhs(ir,2) + flux(2);
   end
end %e