%----------------------------------------------------------------------%
%This subroutine builds the volume integral contribution using the Weak Form DGM-SEM
%on Quadrilateral Elements.
%Written by Francis X. Giraldo on 1/2001
%           Naval Postgraduate School
%           Department of Applied Mathematics
%           Monterey, CA 93943-5502
%----------------------------------------------------------------------%
function rhs = create_rhs_volume(q,u,v,ksi_x,ksi_y,eta_x,eta_y,jac,...
		 wnq,psi,dpsi,intma,iperiodic,npoin,nelem,ngl,form_method)

%global arrays
rhs=zeros(npoin,1);

%Construct Volume Integral Contribution
for e=1:nelem

    %Loop Integration Points
    for j=1:ngl
        for i=1:ngl
            wq=wnq(i)*wnq(j)*jac(i,j,e);
            e_x=ksi_x(i,j,e);
            e_y=ksi_y(i,j,e);
            n_x=eta_x(i,j,e);
            n_y=eta_y(i,j,e);

            %Interpolate at Integration Point
            I=intma(i,j,e);
            u_k=u(I);
            v_k=v(I);
            q_k=q(I);
            f_k=q_k*u_k;
            g_k=q_k*v_k;

            if strcmp(form_method,'weak')
                %Interpolate at Integration Points
                for k=1:ngl
                    dhqde=dpsi(k,i)*psi(j,j)*( e_x*f_k + e_y*g_k );
                    I=iperiodic(intma(k,j,e));
                    rhs(I)=rhs(I) + wq*dhqde;
                    
                    dhqdn=psi(i,i)*dpsi(k,j)*( n_x*f_k + n_y*g_k );
                    I=iperiodic(intma(i,k,e));
                    rhs(I)=rhs(I) + wq*dhqdn;
                end %k
            elseif strcmp(form_method,'strong')
                %Compute derivs in Computational Space
                f_e=0; f_n=0; g_e=0; g_n=0;
                for k=1:ngl
                    h_e=dpsi(k,i)*psi(j,j);
                    I=intma(k,j,e);
                    f_e=f_e + h_e*u(I)*q(I);
                    g_e=g_e + h_e*v(I)*q(I);
                    h_n=psi(i,i)*dpsi(k,j);
                    I=intma(i,k,e);
                    f_n=f_n + h_n*u(I)*q(I);
                    g_n=g_n + h_n*v(I)*q(I);
                end %k

                %Compute Derivs in Physical Space
                dfdx_k=f_e*e_x + f_n*n_x;
                dgdy_k=g_e*e_y + g_n*n_y;
                
                %Loop through I points
                I=iperiodic(intma(i,j,e));
                rhs(I)=rhs(I) - wq*(dfdx_k + dgdy_k);

            end %form_method
        end %i
    end %j
end %e