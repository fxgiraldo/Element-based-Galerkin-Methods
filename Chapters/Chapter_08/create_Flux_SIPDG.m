%---------------------------------------------------------------------%
%This function computes the Laplacian Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_Flux_SIPDG(rhs,intma,coord,nelem,ngl,nq,psi,dpsi,q,qe,qe_x,tau,sipdg_flag)

%NIPDG Flag
ipdg_flag=+1; %=-1 gives NIPDG, =+1 gives SIPDG

%FLUX
for e=1:nelem
    
   %-----------------------I=0 Edge-------------------------%
   %Left Edge
   dx=coord(ngl,e)-coord(1,e);
   jac=dx/2;
   ksi_x=2/dx;
   IL=intma(1,e);
   dhdx_i=dpsi(1,1)*ksi_x;
   h_i=psi(1,1);
   q_L=q(IL);
   dqdx_L=0;
   for i=1:ngl
       ip=intma(i,e);
       dqdx_L=dqdx_L + dpsi(i,1)*q(ip)*ksi_x;
   end
   n_L=-1; %(pointing from Left -> Right)
   
   %Right Edge
   if (e>1)
       dx=coord(ngl,e-1)-coord(1,e-1);
       ksi_x=2/dx;
       IR=intma(ngl,e-1);
       q_R=q(IR);
       dqdx_R=0;
       for i=1:ngl
           ip=intma(i,e-1);
           dqdx_R=dqdx_R + dpsi(i,nq)*q(ip)*ksi_x;
       end
   else %Left BC
       q_R=qe(intma(1,1));
       dqdx_R=qe_x(intma(1,1));
%        q_L=q_R;
%        dqdx_L=dqdx_R;
   end
   q_star=0.5*(q_L + q_R);
   mu=(ngl-1)*(ngl-1+1)/(dx)*tau; %Shabazi
   dqdx_star=0.5*( dqdx_L + dqdx_R)*sipdg_flag - mu*(q_R-q_L);
   
   %Add to Element
   rhs(IL)=rhs(IL) + n_L*h_i*dqdx_star + ipdg_flag*n_L*dhdx_i*(q_L-q_star)*sipdg_flag;

   %-----------------------I=N Edge-------------------------%
   %Left Edge
   IL=intma(ngl,e);
   dhdx_i=dpsi(ngl,nq)*ksi_x;
   h_i=psi(ngl,nq);
   q_L=q(IL);
   dqdx_L=0;
   for i=1:ngl
       ip=intma(i,e);
       dqdx_L=dqdx_L + dpsi(i,nq)*q(ip)*ksi_x;
   end
   n_L=+1; %(pointing from Left -> Right)
   
   %Right Edge
   if (e<nelem)
       dx=coord(ngl,e+1)-coord(1,e+1);
       ksi_x=2/dx;
       IR=intma(1,e+1);
       q_R=q(IR);
       dqdx_R=0;
       for i=1:ngl
           ip=intma(i,e+1);
           dqdx_R=dqdx_R + dpsi(i,1)*q(ip)*ksi_x;
       end
   else %Right BC
       q_R=qe(intma(ngl,nelem));
       dqdx_R=qe_x(intma(ngl,nelem));
%        q_L=q_R;
%        dqdx_L=dqdx_R;
   end
   q_star=0.5*(q_L + q_R);
   mu=(ngl-1)*(ngl-1+1)/(dx)*tau; %Shabazi
   dqdx_star=0.5*( dqdx_L + dqdx_R)*sipdg_flag - mu*(q_R-q_L);
   
   %Add to Element
   rhs(IL)=rhs(IL) + n_L*h_i*dqdx_star + ipdg_flag*n_L*dhdx_i*(q_L-q_star)*sipdg_flag;
end