%----------------------------------------------------------------------%
%This subroutine builds the Normals for a 
%CGDG Method on  Quads
%Written by Francis X. Giraldo
%           Naval Postgraduate School
%           Department of Applied Mathematics
%           Monterey, CA 93943-5216
%----------------------------------------------------------------------%
function [normals,jac_face]=compute_normals(face,intma,coord,...
                          nside,ngl,nq,psi,dpsi)
%global arrays
normals=zeros(2,nq,nside); 
jac_face=zeros(nq,nside);

%local arrays
x=zeros(ngl,ngl);
y=zeros(ngl,ngl);

%loop thru the sides
for is=1:nside

    %Store Left and Right Elements
    ilocl=face(1,is);
    ilocr=face(2,is);
    iel=face(3,is);
    ier=face(4,is);

    %Store Element Variables
    for j=1:ngl
    for i=1:ngl
        ip=intma(i,j,iel);
        x(i,j)=coord(1,ip);
        y(i,j)=coord(2,ip);
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
               normals(1,l,is)=+y_ksi(i,j);
               normals(2,l,is)=-x_ksi(i,j);
               jac_face(l,is)=sqrt( x_ksi(i,j)^2 + y_ksi(i,j)^2 );
            case (2)
               %Side 2: ksi=+1
               i=nq;
               j=l;
               normals(1,l,is)=+y_eta(i,j);
               normals(2,l,is)=-x_eta(i,j);
               jac_face(l,is)=sqrt( x_eta(i,j)^2 + y_eta(i,j)^2 );
            case (3)
               %Side 3: eta=+1
               i=nq+1-l;
               j=nq;
               normals(1,l,is)=-y_ksi(i,j);
               normals(2,l,is)=+x_ksi(i,j);
               jac_face(l,is)=sqrt( x_ksi(i,j)^2 + y_ksi(i,j)^2 ); 
            case (4)
               %Side 4: ksi=-1
               i=1;
               j=nq+1-l;
               normals(1,l,is)=-y_eta(i,j);
               normals(2,l,is)=+x_eta(i,j);
               jac_face(l,is)=sqrt( x_eta(i,j)^2 + y_eta(i,j)^2 );
         end %switch 
    end %l   

    %Normalize Normals
    for l=1:nq
        rnx=sqrt( normals(1,l,is)^2 + normals(2,l,is)^2 );
        normals(1,l,is)=normals(1,l,is)/rnx;
        normals(2,l,is)=normals(2,l,is)/rnx;
     end %l   

end %is
