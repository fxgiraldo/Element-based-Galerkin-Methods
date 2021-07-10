%---------------------------------------------------------------------%
%This function computes the Initial and Exact Solutions.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qe,qb,gravity] = exact_solution(intma,coord,npoin,nelem,ngl,time,icase,eps)

%Initialize
qe=zeros(npoin,2);
qb=zeros(npoin,1);

%Set some constants
xmin=min(min(coord));
xmax=max(max(coord));
xm=0.5*(xmax + xmin);
xc=xm;
xl=xmax-xmin;
amp=0.1;
strength=8;

%Generate Grid Points
for e=1:nelem
   for i=1:ngl
      x=coord(i,e);
      ip=intma(i,e);
      r=x-xc;
      
      if (icase == 1) %Gaussian IC with flat bottom
       gravity=10;
       hmean=0.1;
       hb=0;
       amp=0.5;
       hh=amp*exp( -strength*r^2 );
       qe(ip,1)=hh;
       qe(ip,2)=0;   
       qb(ip)=hmean;
      elseif (icase == 2) %Gaussian IC with Linear Bottom
       gravity=10;
       hmean=0.2;
       hh=amp*exp( -strength*r^2 );
       hb=hmean - 0.1*(1 - (x-xmin)/(xmax-xmin));
       qe(ip,1)=hh;
       qe(ip,2)=0;   
       qb(ip)=hb;
      elseif (icase == 3) %Gaussian IC with Parabolic Bottom
       gravity=10;
       hmean=0.2;
       hh=amp*exp( -strength*r^2 );
       hb=hmean - max(0, 0.1 - 5*x^2);
       qe(ip,1)=hh;
       qe(ip,2)=0;   
       qb(ip)=hb;
      elseif (icase == 4) %Standing Wave with Analytic Solution
       c=2;
       gravity=1;
       hmean=1.0;
       hh=0.5*cos(c*pi*x)*cos(c*pi*time);
       uu=0.5*sin(c*pi*x)*sin(c*pi*time);
       hb=hmean;
       qe(ip,1)=hh;
       qe(ip,2)=uu;   
       qb(ip)=hb;
       elseif (icase == 5) %Dam-Break Problem
          gravity=1;
          hb=0;
          if x <= xc
              hh=1;
              uu=0;
          elseif x > xc
              hh=0.125;
              uu=0;
          end         
          qe(ip,1)=hh;
          qe(ip,2)=uu;   
          qb(ip)=hb;
       elseif (icase == 6) %Dam-Break Problem
          gravity=9.81;
          hb=0;
          cm=1.848; %Simone solution
          cm=4.2585; %FXG solution: between a_l and a_r
          h_l=3;
%           h_r=0.025;
          h_r=1;
          a_l=sqrt(gravity*h_l);
          a_r=sqrt(gravity*h_r);
          xa=xm - time*a_l;
          xb=xm + time*( 2*a_l - 3*cm );
          xc=xm + time*(2*cm^2)*( a_l - cm )/(cm^2 - a_r^2);
          
          a=gravity*h_r;
          b=a_l;
          a6=1;
          a5=0;
          a4=-9*a;
          a3=16*a*b;
          a2=-(8*a*b^2 + a^2);
          a1=0;
          a0=a^3;
          p=[a6 a5 a4 a3 a2 a1 a0];
          r=roots(p);
          
          if x <= xa
              hh=h_l;
              uu=0;
          elseif (x >= xa && x <= xb)
              hh=4/(9*gravity)*( a_l - (x-xm)/(2*time) )^2;
              uu=2/3*( (x-xm)/time + a_l );
          elseif (x >= xb && x <= xc)
              hh=cm^2/gravity;
              uu=2*(a_l - cm);
          elseif (x >= xc)
              hh=h_r;
              uu=0;
          end         
          qe(ip,1)=hh;
          qe(ip,2)=hh*uu;   
          qb(ip)=hb;
      elseif (icase == 7) %Rupert Klein's Linear Multiscale problem
           gravity=1;
           hb=0;
           x0=0.75;
           sigma0=0.1;
           x1=0.25;
           k=7*2*pi;
           p0=exp(-((x-x0)/sigma0)^2);
           p1=exp(-((x-x1)/sigma0)^2)*cos( k*(x-x1)/sigma0 );
           hh=p0 + p1;
           qe(ip,1)=hh;
           qe(ip,2)=hh; 
           qb(ip)=hb;
      end
   end %i
end %e  


      
