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
      
      if (icase == 1) %Gaussian IC with flat bottom
       gravity=10;
       hmean=0.1;
       hb=0;
       amp=0.5;
       hh=amp*exp( -strength*r^2 );
       qe(1,i,ie)=hh;
       qe(2,i,ie)=0;   
       qb(i,ie)=hmean;
      elseif (icase == 2) %Gaussian IC with Linear Bottom
       gravity=10;
       hmean=0.2;
       hh=amp*exp( -strength*r^2 );
       hb=hmean - 0.1*(1 - (x-xmin)/(xmax-xmin));
       qe(1,i,ie)=hh;
       qe(2,i,ie)=0;   
       qb(i,ie)=hb;
      elseif (icase == 3) %Gaussian IC with Parabolic Bottom
       gravity=10;
       hmean=0.2;
       hh=amp*exp( -strength*r^2 );
       hb=hmean - max(0, 0.1 - 5*x^2);
       qe(1,i,ie)=hh;
       qe(2,i,ie)=0;   
       qb(i,ie)=hb;
      elseif (icase == 4) %Standing Wave with Analytic Solution
       c=2;
       gravity=1;
       hmean=1.0;
       hh=0.5*cos(c*pi*x)*cos(c*pi*time);
       uu=0.5*sin(c*pi*x)*sin(c*pi*time);
       hb=hmean;
       qe(1,i,ie)=hh;
       qe(2,i,ie)=uu;   
       qb(i,ie)=hb;
      end
   end %i
end %ie  


      
