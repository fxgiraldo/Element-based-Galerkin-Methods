%---------------------------------------------------------------------%
%This function computes the Initial and Analytic Solutions.
%Written by F.X. Giraldo on 5/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qe] = exact_solution(coord,npoin,time,icase,eq_set)

%Initialize
qe=zeros(3,npoin);

%Set some constants
qb=1.0;
w=1;
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

%Compute Initial Condition
for I=1:npoin
   x=coord(1,I);
   y=coord(2,I);
   r=sqrt( (x-xc)^2 + (y-yc)^2 );
   xx=x - xc*cos(time) - yc*sin(time);
   yy=y + xc*sin(time) - yc*cos(time);
   if strcmp(eq_set,'advection') 
        if (icase == 1)  %Gaussian in CCW direction
          qe(1,I)=a*exp( - sigma*(xx^2 + yy^2 ) );
          qe(2,I)=+w*(y-ym);
          qe(3,I)=-w*(x-xm);
        elseif (icase == 2) %Gaussian along X
          qe(1,I)=a*exp( - sigma*(x^2 + y^2 ) );
          qe(2,I)=w*xl;
          qe(3,I)=0;
        elseif (icase == 3) %Gaussian along Y
          qe(1,I)=a*exp( - sigma*(x^2 + y^2 ) );
          qe(2,I)=0;
          qe(3,I)=w*yl;
        elseif (icase == 4) %Gaussian along X-Y Diagonal
          qe(1,I)=a*exp( - sigma*(x^2 + y^2 ) );
          qe(2,I)=w*xl;
          qe(3,I)=w*yl;
        elseif (icase == 5) %Square in CCW direction
          qe(1,I)=0;
          if (r <= rc)
             qe(1,I)=1.0;
          end
          qe(2,I)=+w*(y-ym);
          qe(3,I)=-w*(x-xm);
        elseif (icase == 6) %Square along X
          qe(1,I)=0;
          if (r <= rc)
             qe(1,I)=1.0;
          end
          qe(2,I)=w*xl;
          qe(3,I)=0;
        else  
          disp(['This IC is not supported: eq_set= ',num2str(eq_set),' icase= ',num2str(icase) ]);
          disp(['Program paused: kill it with ctrl c' ]);
          pause(30);
        end
   elseif strcmp(eq_set,'shallow') 
       if (icase == 7) %Riemann problem for SWE
          qe(1,I)=a*exp( - sigma*(x^2 + y^2 ) );
          qe(2,I)=0;
          qe(3,I)=0;
       else  
          disp(['This IC is not supported: eq_set= ',num2str(eq_set),' icase= ',num2str(icase) ]);
          disp(['Program paused: kill it with ctrl c' ]);
          pause(30);
       end 
   end

   %Compute conservation variables
   qe(1,I)=qe(1,I) + qb;
   qe(2,I)=qe(1,I)*qe(2,I);
   qe(3,I)=qe(1,I)*qe(3,I);
end %I    



      
