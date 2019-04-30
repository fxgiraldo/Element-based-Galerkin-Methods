%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Naval Research Laboratory 
%           Monterey, CA 93943-5502
%---------------------------------------------------------------------%
function fe = source_function(coord,npoin,time,icase)

%Set some constants
w=1;
visc=0;
h=0;
xmin=-1;
xmax=+1;
xl=xmax-xmin;
sigma0=0.125;
xc=-0.7;
%xc=-0.75;
rc=1.e-6;
%rc=0.125;
fc=10;
sigma=sqrt( sigma0^2 + 2*visc*time );
u=w*xl;
icase;

%Initialize
fe=zeros(npoin,1);

timec=time - floor(time);

%Generate Grid Points
for i=1:npoin
  x=coord(i);
  r=abs(x-xc);
  if (icase <= 2) 
     fe(i)=0;
  elseif (icase == 3)
     if (r <= rc) 
       fe(i)=fc;
     end
  elseif (icase == 4)
     if (r <= rc) 
       fe(i)=fc;
     end
  elseif (icase == 5) 
     fe(i)=0;
  end
end %i     


      
