%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Maths
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_flux(rhs,qp,qb,nelem,ngl,diss,h_eps,gravity,delta_nl,icase,time)

%Integrate Flux Terms
for ie=0:nelem
   iel=ie;
   ier=ie+1;
   
   %------Left Variables
   if (ie > 0)
       ps_l=qp(1,ngl,iel);
       pu_l=qp(2,ngl,iel);
       b_l =qb(ngl,iel);
       p_l=ps_l + b_l;
       u_l =pu_l/(p_l);
   end
    
   %------Right Variables
   if (ie < nelem)
       ps_r=qp(1,1,ier);
       pu_r=qp(2,1,ier);
       b_r =qb(1,ier);
       p_r=ps_r + b_r;
       u_r =pu_r/(p_r);
   end
   
   %------Left NFBC
   if (ie == 0)
        p_l=p_r;
        ps_l=ps_r;
        b_l=b_r;
        u_l=-u_r;
   end
   %------Right NFBC
   if (ie == nelem)
        p_r=p_l;
        ps_r=ps_l;
        b_r=b_l;
        u_r=-u_l;
        
        if (icase == 13 || icase == 14 || icase == 15)
            p_r=p_l;
            b_r=b_l;
            u_r=u_l;
%             ps_r=((2*sin(8*atan(1.0)*time/43200))- 1e-2);
%             ps_r=((2*sin(2*pi*time/7200)) + 1e-2);
            ps_r=((2*sin(2*pi*time/7200)) + 0e-2);
            p_r=ps_r + b_r;
        end
        if (icase == 16)
            p_r=p_l;
            b_r=b_l;
            u_r=u_l;
            ps_r=((0.75*sin(8*atan(1.0)*time/3600)) - 1e-2);
            p_r=ps_r + b_r;
        end
   end

   %------Check for small water height
   if (p_l <= h_eps)
       p_l=h_eps;
       %ps_l=p_l - b_l;
       u_l=0;
   end
   
   if (p_r <= h_eps)
       p_r=h_eps;
       %ps_r=p_r - b_l;
       u_r=0;
   end
   
   %Form big U
   pu_l=p_l*u_l;
   pu_r=p_r*u_r;
   
   %Compute Wave Speeds 
   clam_l=( abs(u_l) + sqrt(gravity*p_l) )*delta_nl + sqrt(gravity*b_l)*(1-delta_nl);
   clam_r=( abs(u_r) + sqrt(gravity*p_r) )*delta_nl + sqrt(gravity*b_r)*(1-delta_nl);
   clam=max(clam_l,clam_r);
   
   %Mass Flux
   f_l=pu_l;
   f_r=pu_r;
   flux_p=0.5*( f_l + f_r - diss*abs(clam)*(ps_r - ps_l) );

   %Momentum Flux
   f_l=( pu_l^2/(p_l) + 0.5*gravity*ps_l^2 )*delta_nl + gravity*ps_l*b_l;
   f_r=( pu_r^2/(p_r) + 0.5*gravity*ps_r^2 )*delta_nl + gravity*ps_r*b_r;
   flux_pu=0.5*( f_l + f_r - diss*abs(clam)*(pu_r - pu_l));

   %Make Sure No flux is going to dry element
   if (flux_p > 0)
        if (((ps_l + b_l) <= h_eps)  )
            flux_p=0;
        end
    end
    
    if (flux_p < 0)
        if (((ps_r + b_r) <= h_eps)  )
            flux_p=0;
        end
    end
    
   %Add RHS to Left
   if (ie > 0)
       rhs(1,ngl,iel)=rhs(1,ngl,iel) - flux_p;
       rhs(2,ngl,iel)=rhs(2,ngl,iel) - flux_pu;
   end
   %Add RHS to Right
   if (ie < nelem)
       rhs(1,1,ier)=rhs(1,1,ier) + flux_p;
       rhs(2,1,ier)=rhs(2,1,ier) + flux_pu;
   end
end %ie

