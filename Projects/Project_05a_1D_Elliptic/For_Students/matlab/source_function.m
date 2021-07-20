%---------------------------------------------------------------------%
%This function computes the RHS source function
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function fe = source_function(coord,npoin,time,icase)

%Set some constants
w=1;
visc=0;
h=0;
xmin=-2;
xmax=+2;
xl=xmax-xmin;
sigma0=0.04;
xc=-1.8;
%rc=1e-6;
rc=0.04;
fc=12;
sigma=sqrt( sigma0^2 + 2*visc*time );
u=1;

%Initialize
fe=zeros(npoin,1);

timec=time - floor(time);

%Generate Grid Points
for ip=1:npoin
   x=coord(ip);
   r=abs(x-xc);
   if (icase <= 2) 
      fe(ip)=0;
   elseif (icase == 3)
      fe(ip)=fc*sigma0/sigma*exp( -(x-xc)^2/(sigma^2) );
      if (r <= rc) 
%        fe(ip)=fc;
      end
   elseif (icase == 4)
      fe(ip)=fc*sigma0/sigma*exp( -(x-xc)^2/(sigma^2) );
      if (r <= rc) 
%        fe(ip)=fc;
      end
   end
end %ip      


      
