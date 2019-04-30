%---------------------------------------------------------------------%
%This function computes the Initial and Exact Solutions.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function qe = exact_solution(coord,npoin,time,icase)

%Set some constants
w=1;
visc=0;
h=0;
xc=0;
xmin=-1;
xmax=+1;
xl=xmax-xmin;
sigma0=0.125;
rc=0.125;
sigma=sqrt( sigma0^2 + 2*visc*time );
u=w*xl;
icase;

%Initialize
qe=zeros(npoin,1);

timec=time - floor(time);

%Generate Grid Points
for i=1:npoin
  x=coord(i);
  xbar=xc + u*timec;
  if (xbar > xmax) 
     xbar=xmin + (xbar-xmax);
  end	 
  r=x-xbar;
  if (icase == 1)
      qe(i)=sigma0/sigma*exp( -(x-xbar)^2/(4*sigma^2) );
%            qe(i,ie)=sigma0/sigma*exp( -(x-xbar)^2/(sigma^2) );
  elseif (icase == 2)
     if ( abs(r) <= rc )
        qe(i)=1;
     end
  elseif (icase == 3)
     qe(i)=sigma0/sigma*exp( -(x-xbar)^2/(2*sigma^2) );
  elseif (icase == 4)
     if ( abs(r) <= rc )
        qe(i)=1;
     end
  elseif (icase == 5)
     if ( x <= xc )
        qe(i)=1;
     end
  end
end %i      


      
