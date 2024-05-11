%----------------------------------------------------------------------%
%This subroutine builds the volume integral contribution using the Weak Form DGM-SEM
%on Quadrilateral Elements.
%Written by Francis X. Giraldo on 7/2021
%           Naval Postgraduate School
%           Department of Applied Mathematics
%           Monterey, CA 93943-5502
%----------------------------------------------------------------------%
function rhs = create_rhs_volume_ES(q,u,v,ksi_x,ksi_y,eta_x,eta_y,jac,...
		 wnq,psi,dpsi,intma,iperiodic,npoin,nelem,ngl,flux_method)

%global arrays
rhs=zeros(npoin,1);
flux_es=zeros(ngl,ngl,ngl,ngl,2);
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

    %Compute Flux Matrix
    for j=1:ngl
        for i=1:ngl
            I=node(i,j);
            q_I=q(I); u_I=u(I); v_I=v(I);
            for l=1:ngl
                for k=1:ngl
                    J=node(k,l);
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
                    flux_es(i,j,k,l,1)=2*q_ave*u_ave;
                    flux_es(i,j,k,l,2)=2*q_ave*v_ave;
                end
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
                %Ksi Terms
                h_e=dpsi(k,i)*psi(j,j);
                f_e=flux_es(i,j,k,j,1);
                g_e=flux_es(i,j,k,j,2);
                flux_e=h_e*( f_e*e_x + g_e*e_y);
                %Eta Terms
                h_n=psi(i,i)*dpsi(k,j);
                f_n=flux_es(i,j,i,k,1);
                g_n=flux_es(i,j,i,k,2);
                flux_n=h_n*( f_n*n_x + g_n*n_y);
                %RHS Contribution
                rhs(I)=rhs(I) - wq*(flux_e + flux_n);                
            end %k
        end %i
    end %j
end %e