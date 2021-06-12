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
xmin=min(coord(1,:)); 
xmax=max(coord(1,:));
ymin=min(coord(2,:)); 
ymax=max(coord(2,:));
xl=xmax-xmin;
yl=ymax-ymin;
xm=0.5*(xmax+xmin);
ym=0.5*(ymax+ymin);
xc=xmin + 0.25*xl;
yc=ymin + 0.5*yl;
rc=0.25;
a=1;
sigma=32.0;
timec=time - floor(time);

%Generate Grid Points
for ip=1:npoin
   x=coord(1,ip);
   y=coord(2,ip);
   r=sqrt( (x-xc)^2 + (y-yc)^2 );
   %if (r <= rc)
      xx=x - xc*cos(time) - yc*sin(time);
      yy=y + xc*sin(time) - yc*cos(time);
      qe(ip)=a*exp( - sigma*(xx^2 + yy^2 ) );
   %end
   if (icase == 1)  %Gaussian in CCW direction
      ue(ip)=+w*(y-ym);
      ve(ip)=-w*(x-xm);
   elseif (icase == 2) %Gaussian along X
      xx=x - 0*xc;
      yy=y - 0*yc;
      ue(ip)=w*xl;
      ve(ip)=0;
      qe(ip)=a*exp( - sigma*(xx^2 + yy^2 ) );
   elseif (icase == 3) %Gaussian along Y
      ue(ip)=0;
      ve(ip)=w*yl;
      qe(ip)=a*exp( - sigma*(x^2 + y^2 ) );
   elseif (icase == 4) %Gaussian along X-Y Diagonal
      ue(ip)=w*xl;
      ve(ip)=w*yl;
      qe(ip)=a*exp( - sigma*(x^2 + y^2 ) );
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



      
