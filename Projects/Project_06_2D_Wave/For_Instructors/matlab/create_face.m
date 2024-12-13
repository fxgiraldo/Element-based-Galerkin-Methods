%----------------------------------------------------------------------%
%This subroutine constructs the Side Information for a High Order 
%CG/DG on Quads
%Written by Francis X. Giraldo 
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%----------------------------------------------------------------------%
function [face,mapL,mapR]=create_face(iside,intma,nface,ngl)

%global arrays
face=zeros(4,nface);
mapL=zeros(2,ngl,4);
mapR=zeros(2,ngl,4);

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
   mapL(1,l,1)=l;
   mapL(2,l,1)=1;
   mapR(1,l,1)=ngl+1-l; %because side 1 on the right is = to side 3 on the left
   mapR(2,l,1)=1;

   %ksi=+1
   mapL(1,l,2)=ngl;
   mapL(2,l,2)=l;
   mapR(1,l,2)=ngl;
   mapR(2,l,2)=ngl+1-l; %because side 2 on the right is = to side 4 on the left

   %eta=+1
   mapL(1,l,3)=ngl+1-l; %because side 3 on the left is = to side 1 on the right
   mapL(2,l,3)=ngl;
   mapR(1,l,3)=l;
   mapR(2,l,3)=ngl;

   %ksi=-1
   mapL(1,l,4)=1;
   mapL(2,l,4)=ngl+1-l; %because side 4 on the left is = to side 2 on the right
   mapR(1,l,4)=1;
   mapR(2,l,4)=l;
end %l

%loop thru the sides
for i=1:nface

   ip1=iside(1,i);
   ip2=iside(2,i);
   iel=iside(3,i);
   ier=iside(4,i);

   %check for position on Left Element
   for j=1:4
      j1=j;
      j2=j+1;
      if (j2 > 4) 
         j2=1;
      end %j2   
      jp1=intma(inode(j1),jnode(j1),iel);
      jp2=intma(inode(j2),jnode(j2),iel);

      if (ip1 == jp1 && ip2 == jp2) 
         face(1,i)=j;
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
         jp1=intma(inode(j1),jnode(j1),ier);
         jp2=intma(inode(j2),jnode(j2),ier);

         if (ip1 == jp2 && ip2 == jp1) 
            face(2,i)=j;
            break;          %leave J loop
         end %ip1
      end %j
   end %ier
   
   %Store Elements into FACE
   face(3,i)=iel;
   face(4,i)=ier;
end %i



