%---------------------------------------------------------------------%
%This function computes the Initial and Exact Solutions.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qe,qb,gravity] = exact_solution_dg(coord,nelem,ngl,time,icase,eps)

%Set some constants
xc=0;
xmin=min(min(coord));
xmax=max(max(coord));
xl=xmax-xmin;
amp=0.1;
strength=8;

%Initialize
qe=zeros(2,ngl,nelem);
qb=zeros(ngl,nelem);

%Generate Grid Points
for ie=1:nelem
    for i=1:ngl
        x=coord(i,ie);
        r=x-xc;
        
        if (icase == 1) %Gaussian wave with flat bathymetry
            gravity=10;
            hmean=0.01;
            hb=0;
            amp=0.5;
            hh=amp*exp( -strength*r^2 );
            qe(1,i,ie)=hh;
            qe(2,i,ie)=0;
            qb(i,ie)=1;
        elseif (icase == 20) %Still water  with flat bathymetry
            gravity=10;
            hmean=1;
            qe(1,i,ie) = hmean;
            qe(2,i,ie) = 0;
            qb(i,ie)=0;
        elseif (icase == 2) %Gaussian Wave with linear bathymetry
            gravity=10;
            hmean=0.2;
            %amp=0;
            hh=amp*exp( -strength*r^2 );
            hb=hmean - 0.1*(1 - (x-xmin)/(xmax-xmin));
            qe(1,i,ie)=(hh+0*hmean);
            qe(2,i,ie)=0;
            qb(i,ie)=hb;
        elseif (icase == 3) %Gaussian wave with Quadratic bathymetry
            gravity=10;
            hmean=0.2;
            %amp=0;
            hh=amp*exp( -strength*r^2 );
            hb=hmean - max(0, 0.1 - 5*x^2);
            %hb=max(0, 0.25 - 5*x^2);
            qe(1,i,ie)=(hh+0*hmean);
            qe(2,i,ie)=0;
            qb(i,ie)=hb;
        elseif (icase == 4) %Gaussian wave with Linear Bathymetry
            gravity=10;
            amp=0.2;
            hh=amp*exp( -strength*r^2 );
            hb=1*(x-xmin)/(xmax-xmin) + 1;
            hh=max(1.5-hb,hh);
            qe(1,i,ie)=gravity*hh;
            qe(2,i,ie)=0;
            qb(i,ie)=gravity*hb-15;
        elseif (icase ==5) %WD with Quadratic bathymetry
            gravity=10.0;
            hbtmp=-1e-1+(x)^2;
            hb=min(0.125,hbtmp);
            hh=max(0.0,-hb); %still water condition
            %hh=0;
            qe(1,i,ie)=gravity*hh;
            qe(2,i,ie)=0;
            qb(i,ie)=gravity*hb;
        elseif (icase == 6) %WD with Quadratic bathymetry
            gravity=10.0;
            hb=.5-2*(x^2);
            hh=max(x,-hb);
            qe(1,i,ie)=hh;
            qe(2,i,ie)=0;
            qb(i,ie)=hb;
        elseif (icase == 7) %WD with Linear bathymetry
            gravity=10.0;
            hb=x;
            hh=max(0.0,-hb);
            qe(1,i,ie)=hh;
            qe(2,i,ie)=0;
            qb(i,ie)=hb;
        elseif (icase == 8) %WD with Linear Bathymetry
            gravity=10.0;
            hb=-x*10;
            hh=max(0.0,-hb);
            qe(1,i,ie)=hh;
            qe(2,i,ie)=0;
            qb(i,ie)=hb;
        elseif(icase == 9) %WD with Linear Island
            gravity=10.0;
            hb=max(x,-x);
            hh=max(-.25,-hb);
            qe(1,i,ie)=hh;
            qe(2,i,ie)=0;
            qb(i,ie)=hb;
        elseif(icase == 10) %WD with Quadratic Beach: Right Side
            gravity=10.0;
            hb=min((x-1)^2,2);
            hh=max(-1.0,-hb);
            qe(1,i,ie)=hh;
            qe(2,i,ie)=0;
            qb(i,ie)=hb;
        elseif(icase == 11) %WD with Quadratic Beach: Left Side
            gravity=10.0;
            hb=min((x+1)^2,2);
            hh=max(-1.0,-hb);
            qe(1,i,ie)=hh;
            qe(2,i,ie)=0;
            qb(i,ie)=hb;
        elseif (icase == 12) %WD with Quadratic Island: No Dry land
            gravity=10.0;
            hb=x^2;
            hh=max(.25,-hb);
            qe(1,i,ie)=hh;
            qe(2,i,ie)=0;
            qb(i,ie)=hb;
        elseif (icase == 13) %Balzano Test 1
            gravity=10;
%             hh=0.01;
            hh=0.0;
            hb=x/2760;
            qe(1,i,ie)=hh;
            qe(2,i,ie)=0;
            qb(i,ie)=hb;
        elseif (icase == 14) %Balzano Test 2
            gravity=10;
            hh=0.43;
            if(x < 3600)
                hb=x/2760;
            elseif (x < 4800)
                hb=30/23;
            elseif (x < 6000)
                hb=x/1380 - 50/23;
            else
                hb=x/2760;
            end
            qe(1,i,ie)=hh;
            qe(2,i,ie)=0;
            qb(i,ie)=hb;
        elseif (icase == 15) %Balzano Test 3
            gravity=10;
            hh=0.01;
            %hh=0;
            if(x < 3600)
                hb=x/2760;
            elseif (x < 4800)
                hb=-x/2760+ 60/23;
            elseif (x < 6000)
                hb=x/920 - 100/23;
            else
                hb=x/2760;
            end
            qe(1,i,ie)=hh;
            qe(2,i,ie)=0;
            qb(i,ie)=hb;
        elseif (icase == 16) %WD with Linear Bathymetry
            gravity=10;
            hh=0.01;
            hb=x/1000 + 0.5;
            if( x < 100)
                hb=x/1000 - 0.4;
            elseif(x < 200)
                hb = x/100 - 1.3;
            end
            qe(1,i,ie)=hh;
            qe(2,i,ie)=0;
            qb(i,ie)=hb;
        end
    end %i
end %ie

%Reset Total Water Height
%qe(1,:,:)=qe(1,:,:)-qb(:,:);


