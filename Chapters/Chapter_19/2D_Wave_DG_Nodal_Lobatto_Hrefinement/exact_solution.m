%---------------------------------------------------------------------%
%This function computes the Initial and Analytic Solutions.
%Written by F.X. Giraldo on 5/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qe,ue,ve] = exact_solution(coord,npoin,time,icase)

%Initialize
qe=zeros(npoin,1);
ue=zeros(npoin,1);
ve=zeros(npoin,1);

%Set some constants
w=1;
visc=0;
xmin=min(coord(:,1)); 
xmax=max(coord(:,1));
ymin=min(coord(:,2)); 
ymax=max(coord(:,2));
xl=xmax-xmin;
yl=ymax-ymin;
xm=0.5*(xmax+xmin);
ym=0.5*(ymax+ymin);
xc=xmin + 0.25*xl;
yc=ymin + 0.5*yl;
sigma0=0.125*0.5;
sigma0=1.0/16.0;
rc=0.25;
sigma=sqrt( sigma0^2 + 2*visc*time );
a=0.5*1/(1 + 2*visc*time/sigma0^2);
den=2*(sigma0^2 + 2*visc*time);
timec=time - floor(time);

%Generate Grid Points
for ip=1:npoin
   x=coord(ip,1);
   y=coord(ip,2);
   r=sqrt( (x-xc)^2 + (y-yc)^2 );
   %if (r <= rc)
      xx=x - xc*cos(time) - yc*sin(time);
      yy=y + xc*sin(time) - yc*cos(time);
      qe(ip)=a*exp( - (xx^2 + yy^2 )/den );
   %end
   if (icase == 1)  %Gaussian in CCW direction
      ue(ip)=+w*(y-ym);
      ve(ip)=-w*(x-xm);
      qe(ip)=exp( - 128*(xx^2 + yy^2 ) );
   elseif (icase == 2) %Gaussian along X
      xx=x - 0*xc;%xx=x-2*time;
%       if(xx<xmin)
%          xx = xx+xmax-xmin; 
%       end
      yy=y - 0*yc;
      ue(ip)=w*xl;
      ve(ip)=0;
      qe(ip)=a*exp( - (xx^2 + yy^2 )/den );
   elseif (icase == 3) %Gaussian along Y
      ue(ip)=0;
      ve(ip)=w*yl;
      qe(ip)=a*exp( - (x^2 + y^2 )/den );
   elseif (icase == 4) %Gaussian along X-Y Diagonal
      ue(ip)=w*xl;
      ve(ip)=w*yl;
      qe(ip)=a*exp( - ((x+0.3)^2 + (y+0.3)^2 )/den );
   elseif (icase == 5) %Square in CCW direction
      qe(ip)=0;
      if (r <= rc)
         qe(ip)=1.0;
      end
      ue(ip)=+w*(y-ym);
      ve(ip)=-w*(x-xm);
   elseif (icase == 6) %Square along X
      qe(ip)=0;
      if (r <= rc)
         qe(ip)=1.0;
      end
      ue(ip)=w*xl;
      ve(ip)=0;
   end
end %ip  



      
