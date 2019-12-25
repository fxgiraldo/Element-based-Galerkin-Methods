%----------------------------------------------------------------------%
%This subroutine constructs the Side Information for a High Order 
%Spectal Element Quads
%Written by Francis X. Giraldo 
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%----------------------------------------------------------------------%
function [psideh,imapl,imapr]=create_side_dg(iside,intma,nside,nelem,ngl,nbound,bsidoa,jeside)

%global arrays
psideh=zeros(nside,6);
imapl=zeros(4,2,ngl);
imapr=zeros(4,2,ngl);

%local arrays
inode=zeros(4,1);
jnode=zeros(4,1);

%Construct Boundary Pointer
inode(1)=1;
inode(2)=ngl;
inode(3)=ngl;
inode(4)=1;
jnode(1)=1;
jnode(2)=1;
jnode(3)=ngl;
jnode(4)=ngl;

%Construct IMAP arrays
for l=1:ngl

   %eta=-1
   imapl(1,1,l)=l;
   imapl(1,2,l)=1;
%    imapr(1,1,l)=ngl+1-l;
   imapr(1,1,l)=l;
   imapr(1,2,l)=1;

   %ksi=+1
   imapl(4,1,l)=ngl;
   imapl(4,2,l)=l;
   imapr(4,1,l)=ngl;
   imapr(4,2,l)=l;
%    imapr(2,2,l)=ngl+1-l;

   %eta=+1
%    imapl(3,1,l)=ngl+1-l;
%    imapl(3,2,l)=ngl;
   imapl(2,1,l)=l;
   imapl(2,2,l)=ngl;
   imapr(2,1,l)=l;
   imapr(2,2,l)=ngl;

   %ksi=-1
   imapl(3,1,l)=1;
%    imapl(4,2,l)=ngl+1-l;
   imapl(3,2,l)=l;
   imapr(3,1,l)=1;
   imapr(3,2,l)=l;
end %l

%loop thru the sides
for i=1:nside

   ip1=iside(i,1);
   ip2=iside(i,2);
   iel=iside(i,3);
   ier=iside(i,4);

   %check for position on Left Element
   for j=1:4
      j1=j;
      j2=j+1;
      if (j2 > 4) 
         j2=1;
      end %j2   
      jp1=intma(iel,inode(j1),jnode(j1));
      jp2=intma(iel,inode(j2),jnode(j2));

      if (ip1 == jp1 && ip2 == jp2) 
         psideh(i,1)=j;
         break; %leave J loop
      end %ip1
   end %j

   %check for position on Right Element
   if (ier > 0)
      for j=1:4
         j1=j;
         j2=j+1;
         if (j2 > 4) 
            j2=1;
         end %j2   
         jp1=intma(ier,inode(j1),jnode(j1));
         jp2=intma(ier,inode(j2),jnode(j2));

         if (ip1 == jp2 && ip2 == jp1) 
            psideh(i,2)=j;
            break;          %leave J loop
         end %ip1
      end %j
   end %ier
   
   if (ier<0)
      for j=1:nbound/2
          if(iel == bsidoa(j,3) && psideh(i,1) == bsidoa(j,1))
             psideh(i,2) = bsidoa(j,2);
             psideh(i,4) = bsidoa(j,4);
             ier = bsidoa(j,4);
             psideh(i,5) = -bsidoa(j,5);
          elseif(iel == bsidoa(j,4) && psideh(i,1) == bsidoa(j,2))
             psideh(i,2) = bsidoa(j,1);
             psideh(i,4) = bsidoa(j,3);
             ier = bsidoa(j,3);
%              psideh(i,5) = -bsidoa(j,5);
          end
          
      end
   end
   
   %Store Elements into PSIDEH
   psideh(i,3)=iel;
   psideh(i,4)=ier;
end %i

%store periodicity information
for i=1:nbound/2
   psideh(jeside(bsidoa(i,3),bsidoa(i,1)),6) = jeside(bsidoa(i,4),bsidoa(i,2));
   psideh(jeside(bsidoa(i,4),bsidoa(i,2)),6) = jeside(bsidoa(i,3),bsidoa(i,1)); 
end



