%---------------------------------------------------------------------%
%This function computes the Courant Number on Quadrilaterals.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [Courant,vel,ds,dt] = compute_Courant(u,v,intma,coord,ngl,dt,Courant_max,sfc,nsfc)

%Initialize
Courant=-1000;
x_min=coord(2,1)-coord(1,1);
y_min=coord(2,2)-coord(1,2);

for is=1:nsfc
      ie=sfc(is); 
    %Loop through I points
    for j=1:ngl-1
    for i=1:ngl-1
        i1=intma(ie,i,j);
        i2=intma(ie,i+1,j);
        i3=intma(ie,i,j+1);
        i4=intma(ie,i+1,j+1);
        x1=coord(i1,1); y1=coord(i1,1);
        x2=coord(i2,1); y2=coord(i2,2);
        x3=coord(i3,1); y3=coord(i3,2);
        x4=coord(i4,1); y4=coord(i4,2);
        ds=sqrt( (x4-x1)^2 + (y4-y1)^2 );
        
        u1=u(ie,i,j); v1=v(ie,i,j);
        u2=u(ie,i+1,j); v2=v(ie,i+1,j);
        u3=u(ie,i,j+1); v3=v(ie,i,j+1);
        u4=u(ie,i+1,j+1); v4=v(ie,i+1,j+1);
        uu=(u1 + u2 + u3 + u4)/4;
        vv=(v1 + v2 + v3 + v4)/4;
        
        vel=sqrt( uu^2 + vv^2 );
        C=vel*dt/ds;
        Courant=max(Courant,C); 
        x_min1=coord(i2,1)-coord(i1,1);
        y_min1=coord(i4,2)-coord(i1,2);    
        if(x_min1<x_min) 
            x_min = x_min1; 
        end
        if(y_min1<y_min) 
            y_min = y_min1;
        end
    end %i
    end %j
end %ie

u_max=max(max(max(u)));
v_max=max(max(max(v)));

ds=sqrt( x_min^2 + y_min^2 );
vel=sqrt( u_max^2 + v_max^2 );
dt=Courant_max*ds/vel;
Courant=Courant_max;  


      
