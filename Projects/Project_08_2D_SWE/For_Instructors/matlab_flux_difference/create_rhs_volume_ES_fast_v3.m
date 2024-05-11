%----------------------------------------------------------------------%
%This subroutine builds the volume integral contribution using the Weak Form DGM-SEM
%on Quadrilateral Elements.
%Written by Francis X. Giraldo on 7/2021
%           Naval Postgraduate School
%           Department of Applied Mathematics
%           Monterey, CA 93943-5502
%----------------------------------------------------------------------%
function rhs = create_rhs_volume_ES_fast_v3(q,u,v,ksi_x,ksi_y,eta_x,eta_y,jac,...
		 wnq,psi,dpsi,intma,iperiodic,npoin,nelem,ngl,flux_method)

%global arrays
rhs=zeros(npoin,1);
flux_e=zeros(ngl,ngl,ngl,2);
flux_n=zeros(ngl,ngl,ngl,2);
node=zeros(ngl,ngl);

%Construct Volume Integral Contribution
for e=1:nelem

    %Loop through element DOF and store Node Pointer and Fluxes
    for j=1:ngl
        for i=1:ngl
            I=intma(i,j,e);
            node(i,j)=I;
        end
    end

    %Compute Flux Matrix: i,j,k,j
    for j=1:ngl
        for i=1:ngl
            I=node(i,j);
            q_I=q(I); u_I=u(I); v_I=v(I);
            for k=1:ngl
                %i,j,k,j=J_e
                J=node(k,j);
                q_J=q(J); u_J=u(J); v_J=v(J);
                if strcmp(flux_method,'OR')
                    q_ave=0.5*q_J;
                    u_ave=u_J;
                    v_ave=v_J;
                elseif strcmp(flux_method,'KG')
                    q_ave=0.5*(q_I + q_J);
                    u_ave=0.5*(u_I + u_J);
                    v_ave=0.5*(v_I + v_J);
                end
                flux_e(i,j,k,1)=2*q_ave*u_ave;
                flux_e(i,j,k,2)=2*q_ave*v_ave;
                %i,j,i,k=J_n
                J=node(i,k);
                q_J=q(J); u_J=u(J); v_J=v(J);
                if strcmp(flux_method,'OR')
                    q_ave=0.5*q_J;
                    u_ave=u_J;
                    v_ave=v_J;
                elseif strcmp(flux_method,'KG')
                    q_ave=0.5*(q_I + q_J);
                    u_ave=0.5*(u_I + u_J);
                    v_ave=0.5*(v_I + v_J);
                end
                flux_n(i,j,k,1)=2*q_ave*u_ave;
                flux_n(i,j,k,2)=2*q_ave*v_ave;
            end
        end
    end

    %Loop Integration Points
    for j=1:ngl
        for i=1:ngl
            I=iperiodic(node(i,j));
            wq=wnq(i)*wnq(j)*jac(i,j,e);
            e_x=ksi_x(i,j,e);
            e_y=ksi_y(i,j,e);
            n_x=eta_x(i,j,e);
            n_y=eta_y(i,j,e);
            
            for k=1:ngl
                %Xi Terms
                h_e=dpsi(k,i)*psi(j,j);
                f_e=flux_e(i,j,k,1);
                g_e=flux_e(i,j,k,2);
                F_flux_e=h_e*( f_e*e_x + g_e*e_y);
                %Eta Terms
                h_n=psi(i,i)*dpsi(k,j);
                f_n=flux_n(i,j,k,1);
                g_n=flux_n(i,j,k,2);
                F_flux_n=h_n*( f_n*n_x + g_n*n_y);
                %RHS Contribution
                rhs(I)=rhs(I) - wq*(F_flux_e + F_flux_n);                
            end %k
        end %i
    end %j
end %e