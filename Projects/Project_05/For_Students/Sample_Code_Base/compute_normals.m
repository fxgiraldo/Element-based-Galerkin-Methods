%----------------------------------------------------------------------%
%This subroutine builds the Normals for a 
%Spectral Element Quads
%Written by Francis X. Giraldo on 1/01 (from Brown University)
%           Naval Postgraduate School
%           Department of Applied Mathematics
%           Monterey, CA 93943-5216
%----------------------------------------------------------------------%
function [nx,ny,jac_face]=compute_normals(face,intma,coord,...
                          nside,ngl,nq,psi,dpsi)
%global arrays
nx=zeros(nside,nq); 
ny=zeros(nside,nq);
jac_face=zeros(nside,nq);

%local arrays
x=zeros(ngl,ngl);
y=zeros(ngl,ngl);

%loop thru the sides
for is=1:nside

    %Store Left and Right Elements
    ilocl=face(is,1);
    ilocr=face(is,2);
    iel=face(is,3);
    ier=face(is,4);

    %Store Element Variables
    for j=1:ngl
    for i=1:ngl
        ip=intma(iel,i,j);
        x(i,j)=coord(ip,1);
        y(i,j)=coord(ip,2);
    end %i  
    end %j  

    %Construct Mapping Derivatives: dx/dksi, dx/deta, dy/dksi,
    %dy/deta
    [x_ksi,x_eta]=map_deriv(psi,dpsi,x,ngl,nq);
    [y_ksi,y_eta]=map_deriv(psi,dpsi,y,ngl,nq);
    
    %Compute Normals 
    for l=1:nq
         switch (ilocl)
            case (1)
               %Side 1: eta=-1
               i=l;
               j=1;
               nx(is,l)=+y_ksi(i,j);
               ny(is,l)=-x_ksi(i,j);
               jac_face(is,l)=sqrt( x_ksi(i,j)^2 + y_ksi(i,j)^2 );
            case (2)
               %Side 2: ksi=+1
               i=nq;
               j=l;
               nx(is,l)=+y_eta(i,j);
               ny(is,l)=-x_eta(i,j);
               jac_face(is,l)=sqrt( x_eta(i,j)^2 + y_eta(i,j)^2 );
            case (3)
               %Side 3: eta=+1
               i=nq+1-l;
               j=nq;
               nx(is,l)=-y_ksi(i,j);
               ny(is,l)=+x_ksi(i,j);
               jac_face(is,l)=sqrt( x_ksi(i,j)^2 + y_ksi(i,j)^2 ); 
            case (4)
               %Side 4: ksi=-1
               i=1;
               j=nq+1-l;
               nx(is,l)=-y_eta(i,j);
               ny(is,l)=+x_eta(i,j);
               jac_face(is,l)=sqrt( x_eta(i,j)^2 + y_eta(i,j)^2 );
         end %switch 
    end %l   

    %Normalize Normals
    for l=1:nq
        rnx=sqrt( nx(is,l)^2 + ny(is,l)^2 );
        nx(is,l)=nx(is,l)/rnx;
        ny(is,l)=ny(is,l)/rnx;
     end %l   

end %is
