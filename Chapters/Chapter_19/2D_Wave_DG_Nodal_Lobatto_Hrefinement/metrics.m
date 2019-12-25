%---------------------------------------------------------------------%
%This function computes the Metric Terms for General 2D Quad Grids.
%Written by F.X. Giraldo on 4/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [ksi_x,ksi_y,eta_x,eta_y,jac] = metrics(coord,intma,psi,dpsi,wnq,nelem,ngl,nq,sfc,nsfc)

%Initialize Global Arrays
ksi_x=zeros(nelem,nq,nq);
ksi_y=zeros(nelem,nq,nq);
eta_x=zeros(nelem,nq,nq);
eta_y=zeros(nelem,nq,nq);
jac=zeros(nelem,nq,nq);

%Initialize Local Arrays
x_ksi=zeros(nq,nq);
x_eta=zeros(nq,nq);
y_ksi=zeros(nq,nq);
y_eta=zeros(nq,nq);
x=zeros(ngl,ngl);
y=zeros(ngl,ngl);

%loop thru the elements
for is=1:nsfc
    ie=sfc(is);

    %Store Element Variables
    for j=1:ngl
    for i=1:ngl
       ip=intma(ie,i,j);
       x(i,j)=coord(ip,1);
       y(i,j)=coord(ip,2);
    end %i  
    end %j

    %Construct Mapping Derivatives: dx/dksi, dx/deta, dy/dksi, dy/deta
    [x_ksi,x_eta]=map_deriv(psi,dpsi,x,ngl,nq);
    [y_ksi,y_eta]=map_deriv(psi,dpsi,y,ngl,nq);

    %Construct Inverse Mapping: dksi/dx, dksi/dy, deta/dx, deta/dy
    for j=1:nq
    for i=1:nq
        xjac=x_ksi(i,j)*y_eta(i,j) - x_eta(i,j)*y_ksi(i,j);
        ksi_x(ie,i,j)=+1.0/xjac*y_eta(i,j);
        ksi_y(ie,i,j)=-1.0/xjac*x_eta(i,j);
        eta_x(ie,i,j)=-1.0/xjac*y_ksi(i,j);
        eta_y(ie,i,j)=+1.0/xjac*x_ksi(i,j);
        jac(ie,i,j)=wnq(i)*wnq(j)*abs(xjac);
    end %i
    end %j
end %ie
