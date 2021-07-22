%---------------------------------------------------------------------%
%This function computes the Metric Terms for General 3D Hexahedral Grids.
%
%This function uses Algs. 12.4 and 12.5 in the textbook, which are the cross-product
%metric terms.
%
%Written by F.X. Giraldo on 7/2021
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [ksi_x,ksi_y,ksi_z,eta_x,eta_y,eta_z,zeta_x,zeta_y,zeta_z,jac] = metrics(coord,intma,psi,dpsi,nelem,ngl,nq)

%Initialize Global Arrays
ksi_x=zeros(nq,nq,nq,nelem);
ksi_y=zeros(nq,nq,nq,nelem);
ksi_z=zeros(nq,nq,nq,nelem);
eta_x=zeros(nq,nq,nq,nelem);
eta_y=zeros(nq,nq,nq,nelem);
eta_z=zeros(nq,nq,nq,nelem);
zeta_x=zeros(nq,nq,nq,nelem);
zeta_y=zeros(nq,nq,nq,nelem);
zeta_z=zeros(nq,nq,nq,nelem);
jac=zeros(nq,nq,nq,nelem);

%Initialize Local Arrays
x=zeros(ngl,ngl,ngl);
y=zeros(ngl,ngl,ngl);
z=zeros(ngl,ngl,ngl);

%loop thru the elements
for e=1:nelem

    %Store Element Variables
    for k=1:ngl
    for j=1:ngl
    for i=1:ngl
       ip=intma(i,j,k,e);
       x(i,j,k)=coord(1,ip);
       y(i,j,k)=coord(2,ip);
       z(i,j,k)=coord(3,ip);
    end %i  
    end %j
    end %k

    %Construct Mapping Derivatives: dx/dksi, dx/deta, dy/dksi, dy/deta
    [x_ksi,x_eta,x_zeta]=map_deriv(psi,dpsi,x,ngl,nq);
    [y_ksi,y_eta,y_zeta]=map_deriv(psi,dpsi,y,ngl,nq);
    [z_ksi,z_eta,z_zeta]=map_deriv(psi,dpsi,z,ngl,nq);

    %Construct Inverse Mapping: dksi/dx, dksi/dy, deta/dx, deta/dy
    for k=1:nq
    for j=1:nq
    for i=1:nq    
        xj = ...
        (x_ksi(i,j,k)*y_eta(i,j,k)*z_zeta(i,j,k) - x_ksi(i,j,k)*y_zeta(i,j,k)*z_eta(i,j,k))  ...
        - (y_ksi(i,j,k)*x_eta(i,j,k)*z_zeta(i,j,k) - y_ksi(i,j,k)*x_zeta(i,j,k)*z_eta(i,j,k))  ...
        + (z_ksi(i,j,k)*x_eta(i,j,k)*y_zeta(i,j,k) - z_ksi(i,j,k)*x_zeta(i,j,k)*y_eta(i,j,k));

        ksi_x(i,j,k,e)=  (y_eta(i,j,k)*z_zeta(i,j,k)-y_zeta(i,j,k)*z_eta(i,j,k))/xj;
        ksi_y(i,j,k,e)= -(x_eta(i,j,k)*z_zeta(i,j,k)-x_zeta(i,j,k)*z_eta(i,j,k))/xj;
        ksi_z(i,j,k,e)=  (x_eta(i,j,k)*y_zeta(i,j,k)-x_zeta(i,j,k)*y_eta(i,j,k))/xj;
        eta_x(i,j,k,e)= -(y_ksi(i,j,k)*z_zeta(i,j,k)-y_zeta(i,j,k)*z_ksi(i,j,k))/xj;
        eta_y(i,j,k,e)=  (x_ksi(i,j,k)*z_zeta(i,j,k)-x_zeta(i,j,k)*z_ksi(i,j,k))/xj;
        eta_z(i,j,k,e)= -(x_ksi(i,j,k)*y_zeta(i,j,k)-x_zeta(i,j,k)*y_ksi(i,j,k))/xj;
        zeta_x(i,j,k,e)= (y_ksi(i,j,k)*z_eta(i,j,k) -y_eta(i,j,k)*z_ksi(i,j,k) )/xj;
        zeta_y(i,j,k,e)=-(x_ksi(i,j,k)*z_eta(i,j,k) -x_eta(i,j,k)*z_ksi(i,j,k) )/xj;
        zeta_z(i,j,k,e)= (x_ksi(i,j,k)*y_eta(i,j,k) -x_eta(i,j,k)*y_ksi(i,j,k) )/xj;    
        jac(i,j,k,e)=abs(xj);    
    end %i
    end %j
    end %k


end %e
