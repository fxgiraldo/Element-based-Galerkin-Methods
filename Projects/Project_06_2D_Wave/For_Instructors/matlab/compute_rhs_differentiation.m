%----------------------------------------------------------------------%
%This subroutine builds the RHS vector for the Strong Form DGM-SEM
%on Quadrilateral Elements for the 2D Euler Equations.
%Written by Francis X. Giraldo on 1/2001
%           Naval Postgraduate School
%           Department of Applied Mathematics
%           Monterey, CA 93943-5502
%----------------------------------------------------------------------%
function rhs = compute_rhs_differentiation(q,u,v,ksi_x,ksi_y,eta_x,eta_y,jac,...
		 wnq,dpsi,intma,iperiodic,npoin,nelem,ngl)

%global arrays
rhs=zeros(npoin,1);

%Construct FEM-type Operators
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

            %Interpolate at Integration Points
            for k=1:ngl
                dhqde=dpsi(k,i)*( e_x*f_k + e_y*g_k );
                I=iperiodic(intma(k,j,e));
                rhs(I)=rhs(I) + wq*dhqde;

                dhqdn=dpsi(k,j)*( n_x*f_k + n_y*g_k );
                I=iperiodic(intma(i,k,e));
                rhs(I)=rhs(I) + wq*dhqdn;
            end %k
        end %i
    end %j
end %e