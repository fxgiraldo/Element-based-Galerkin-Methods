%---------------------------------------------------------------------%
%This function computes the Courant Number on Quadrilaterals.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [Courant,vel,ds,dt] = compute_Courant(u,v,intma,coord,nelem,ngl,Courant_max)

%Initialize
Courant=-1000;
dt=1;

for ie=1:nelem
       
    %Loop through I points
    for j=1:ngl-1
    for i=1:ngl-1
        i1=intma(i,j,ie);
        i2=intma(i+1,j,ie);
        i3=intma(i,j+1,ie);
        i4=intma(i+1,j+1,ie);
        x1=coord(1,i1); y1=coord(2,i1);
        x2=coord(1,i2); y2=coord(2,i2);
        x3=coord(1,i3); y3=coord(2,i3);
        x4=coord(1,i4); y4=coord(2,i4);
        ds=sqrt( (x4-x1)^2 + (y4-y1)^2 );
        
        u1=u(i1); v1=v(i1);
        u2=u(i2); v2=v(i2);
        u3=u(i3); v3=v(i3);
        u4=u(i4); v4=v(i4);
        uu=(u1 + u2 + u3 + u4)/4;
        vv=(v1 + v2 + v3 + v4)/4;
        
        vel=sqrt( uu^2 + vv^2 );
        C=vel*dt/ds;
        Courant=max(Courant,C); 
    end %i
    end %j
end %ie

Courant;

u_max=max(u);
v_max=max(v);
x_min=coord(1,2)-coord(1,1);
y_min=coord(2,2)-coord(2,1);
ds=sqrt( x_min^2 + y_min^2 );
vel=sqrt( u_max^2 + v_max^2 );
dt=Courant_max*ds/vel;
Courant=Courant_max;  


      
