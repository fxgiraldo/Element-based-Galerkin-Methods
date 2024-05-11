%----------------------------------------------------------------------%
%This subroutine builds the volume integral contribution using a united
%CG/DG method on Quadrilateral Elements for the 2D Shallow Water Equations.
%Written by Francis X. Giraldo on 11/2022
%           Naval Postgraduate School
%           Department of Applied Mathematics
%           Monterey, CA 93943-5502
%----------------------------------------------------------------------%
function rhs = create_rhs_volume(q,ksi_x,ksi_y,eta_x,eta_y,jac,...
		 wnq,psi,dpsi,intma,iperiodic,npoin,nelem,ngl,form_method)

%global arrays
rhs=zeros(3,npoin);
flux=zeros(3,2);

%Construct Volume Integral Contribution
for e=1:nelem
    
    %Loop Integration Points
    for j=1:ngl
        for i=1:ngl
            wq=wnq(i)*wnq(j)*jac(i,j,e);
            e_x=ksi_x(i,j,e); e_y=ksi_y(i,j,e); 
            n_x=eta_x(i,j,e); n_y=eta_y(i,j,e);

            %Interpolate at Integration Point
            I=intma(i,j,e);
            r_k=q(1,I);
            u_k=q(2,I)/r_k;
            v_k=q(3,I)/r_k;
            flux(1,1)=r_k*u_k;
            flux(1,2)=r_k*v_k;
            flux(2,1)=r_k*u_k*u_k + 0.5*r_k^2;
            flux(2,2)=r_k*u_k*v_k;
            flux(3,1)=r_k*v_k*u_k;
            flux(3,2)=r_k*v_k*v_k + 0.5*r_k^2;

            if strcmp(form_method,'weak')
                %Interpolate at Integration Points
                for k=1:ngl
                    %Xi derivatives
                    I=iperiodic(intma(k,j,e));
                    h_e=dpsi(k,i)*psi(j,j); 
                    for l=1:3
                        dhde=h_e*( e_x*flux(l,1) + e_y*flux(l,2) );
                        rhs(l,I)=rhs(l,I) + wq*dhde;                    
                    end
                    %Eta derivatives
                    I=iperiodic(intma(i,k,e));
                    h_n=psi(i,i)*dpsi(k,j);
                    for l=1:3
                        dhdn=h_n*( n_x*flux(l,1) + n_y*flux(l,2) );
                        rhs(l,I)=rhs(l,I) + wq*dhdn;                    
                    end
               end %k
            elseif strcmp(form_method,'strong')
                %Xi derivatives
                flux_e=zeros(3,2); 
                for k=1:ngl
                    I=intma(k,j,e);
                    r_k=q(1,I);
                    u_k=q(2,I)/r_k;
                    v_k=q(3,I)/r_k;
                    flux(1,1)=r_k*u_k;
                    flux(1,2)=r_k*v_k;
                    flux(2,1)=r_k*u_k*u_k + 0.5*r_k^2;
                    flux(2,2)=r_k*u_k*v_k;
                    flux(3,1)=r_k*v_k*u_k;
                    flux(3,2)=r_k*v_k*v_k + 0.5*r_k^2;
                    h_e=dpsi(k,i)*psi(j,j);
                    for m=1:2
                        for l=1:3
                            flux_e(l,m)=flux_e(l,m) + h_e*flux(l,m);
                        end %l
                    end %m
                end %k

                %Eta derivatives
                flux_n=zeros(3,2); 
                for k=1:ngl
                    I=intma(i,k,e);
                    r_k=q(1,I);
                    u_k=q(2,I)/r_k;
                    v_k=q(3,I)/r_k;
                    flux(1,1)=r_k*u_k;
                    flux(1,2)=r_k*v_k;
                    flux(2,1)=r_k*u_k*u_k + 0.5*r_k^2;
                    flux(2,2)=r_k*u_k*v_k;
                    flux(3,1)=r_k*v_k*u_k;
                    flux(3,2)=r_k*v_k*v_k + 0.5*r_k^2;
                    h_n=psi(i,i)*dpsi(k,j);
                    for m=1:2
                        for l=1:3
                            flux_n(l,m)=flux_n(l,m) + h_n*flux(l,m);
                        end %l
                    end %m
                end %k

                %x,y derivatives
                flux_x=zeros(3,2); flux_y=zeros(3,2);
                for m=1:2
                    for l=1:3
                        flux_x(l,m)=flux_e(l,m)*e_x + flux_n(l,m)*n_x;
                        flux_y(l,m)=flux_e(l,m)*e_y + flux_n(l,m)*n_y;
                    end %l
                end %m
                
                %Loop through I points
                I=iperiodic(intma(i,j,e));
                for l=1:3
                    rhs(l,I)=rhs(l,I) - wq*(flux_x(l,1) + flux_y(l,2));
                end %l
            end %form_method
        end %i
    end %j
end %e