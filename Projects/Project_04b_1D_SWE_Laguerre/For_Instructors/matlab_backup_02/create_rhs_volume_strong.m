%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 1/2024
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_rhs_volume_strong(qp,qb,intma,jac,npoin,nelem,ngl,wnq,dpsi,gravity,delta_nl,Igamma,icase)

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
        dFhdx_k=0; dFUdx_k=0; dBdx_k=0; 
        for j=1:ngl(e)
            J=intma(j,e);
            dhdx_ji=dpsi(j,i,e)*ksi_x;
            %Form fluxes
            hs_k=qp(J,1);
            hb_k=qb(J);
            h_k=hs_k+hb_k;
            U_k=qp(J,2);
            flux_H=U_k;
            if (icase == -1) %Advection
                flux_H=hs_k*U_k;
            end
            flux_U=( U_k^2/(h_k) + 0.5*gravity*hs_k^2 )*delta_nl + gravity*hb_k*hs_k;
            %Form derivative of fluxes
            dFhdx_k=dFhdx_k + flux_H*dhdx_ji;
            dFUdx_k=dFUdx_k + flux_U*dhdx_ji;
            dBdx_k=dBdx_k   + hb_k*dhdx_ji;
        end

        %Add fluxes
        rhs(I,1)=rhs(I,1) - wq*dFhdx_k;
        rhs(I,2)=rhs(I,2) - wq*dFUdx_k;

        %Add Bathymetry gradient to Momentum
        rhs(I,2)=rhs(I,2) + wq*gravity*hs_k*dBdx_k;

        %Add Rayleigh damping term
        rhs(I,1)=rhs(I,1) - wq*Igamma(I)*hs_k;
        rhs(I,2)=rhs(I,2) - wq*Igamma(I)*U_k;
   end %i
end %e