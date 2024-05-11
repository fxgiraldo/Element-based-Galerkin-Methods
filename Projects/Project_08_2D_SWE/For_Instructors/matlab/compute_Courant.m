%---------------------------------------------------------------------%
%This function computes the Courant Number on Quadrilaterals.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [Courant,vel,ds,dt] = compute_Courant(q,intma,coord,nelem,ngl,Courant_max,eq_set)

%Initialize
C_max=-1000;

flag=1.0;
if strcmp(eq_set,'advection') 
    flag=0.0;
end

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
        
        r1=q(1,i1); u1=q(2,i1)/r1; v1=q(3,i1)/r1;
        r2=q(1,i2); u2=q(2,i2)/r2; v2=q(3,i2)/r2;
        r3=q(1,i3); u3=q(2,i3)/r3; v3=q(3,i3)/r3;
        r4=q(1,i4); u4=q(2,i4)/r4; v4=q(3,i4)/r4;
        rr=(r1 + r2 + r3 + r4)/4;
        uu=(u1 + u2 + u3 + u4)/4;
        vv=(v1 + v2 + v3 + v4)/4;
        
        vel=sqrt( uu^2 + vv^2 ) + flag*sqrt(rr);
        C=vel/ds;
        C_max=max(C_max,C); 
    end %i
    end %j
end %ie

dt=Courant_max/C_max;
Courant=C_max*dt;