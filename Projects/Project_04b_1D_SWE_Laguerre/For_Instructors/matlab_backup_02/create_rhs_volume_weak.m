%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 1/2024
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_rhs_volume_weak(qp,qb,intma,jac,npoin,nelem,ngl,wnq,dpsi,gravity,delta_nl,Igamma,icase)

%Initialize
rhs=zeros(npoin,2);

%Integrate Divergence of Flux (Weak Form)
for e=1:nelem
      
   %Jacobians
   jac_e=jac(e);
   ksi_x=1.0/jac_e;
   
   %LGL Integration
   for i=1:ngl(e)
        wq=wnq(i,e)*jac_e;
        
        %Flux
        I=intma(i,e);
        hs_k=qp(I,1);
        hb_k=qb(I);
        h_k=hs_k+hb_k;
        U_k=qp(I,2);
        flux_H=U_k;
        if (icase == -1) %Advection
            flux_H=hs_k*U_k;
        end
        flux_U=( U_k^2/(h_k) + 0.5*gravity*hs_k^2 )*delta_nl + gravity*hb_k*hs_k;
        
        %Form Derivatives
        dBdx_k=0;
        for j=1:ngl(e)
            J=intma(j,e);
            dhdx_ji=dpsi(j,i,e)*ksi_x;
            dBdx_k=dBdx_k + qb(J)*dhdx_ji;
        end

        %Add Bathymetry gradient to Momentum
        rhs(I,2)=rhs(I,2) + wq*gravity*hs_k*dBdx_k;

        %Add Rayleigh damping term
        rhs(I,1)=rhs(I,1) - wq*Igamma(I)*hs_k;
        rhs(I,2)=rhs(I,2) - wq*Igamma(I)*U_k;
        
        %Form RHS
        for j=1:ngl(e)
            J=intma(j,e);
            dhdx_ji=dpsi(j,i,e)*ksi_x;
            rhs(J,1)=rhs(J,1) + wq*dhdx_ji*flux_H;
            rhs(J,2)=rhs(J,2) + wq*dhdx_ji*flux_U;
        end %j
   end %i
end %e