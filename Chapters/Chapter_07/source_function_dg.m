%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Maths
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function fe = source_function_dg(coord,nelem,ngl,time,icase)

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
fe=zeros(ngl,nelem);

timec=time - floor(time);

%Generate Grid Points
for ie=1:nelem
   for i=1:ngl
      x=coord(i,ie);
      r=abs(x-xc);
      if (icase <= 2) 
         fe(i,ie)=0;
      elseif (icase == 3)
%         fe(i,ie)=sigma0/sigma*exp( -(x-xc)^2/(2*sigma^2) );
         if (r <= rc) 
           fe(i,ie)=fc;
         end
      elseif (icase == 4)
%         fe(i,ie)=sigma0/sigma*exp( -(x-xc)^2/(2*sigma^2) );
         if (r <= rc) 
           fe(i,ie)=fc;
         end
      elseif (icase == 5) 
         fe(i,ie)=0;
      end
   end %i
end %ie      


      
