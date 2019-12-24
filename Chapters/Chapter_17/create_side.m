%----------------------------------------------------------------------%
%This subroutine creates the array ISIDE which stores all of
%the information concerning the sides of all the elements.
%Written by Francis X. Giraldo on 1/01
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5502
%
% INPUT LIST: intma = element connectivity
%             bsido = boundary info (which points are on a boundary, which
%                     element it belongs to and the boundary condition).
%             npoin = number of global points
%             nelem = number of elements
%             nboun = number of boundary faces/edges
%             nside=nface are the number of sides/face/edges in the grid
%             ngl = number of points in one direction in an element
%
% OUTPUT LIST: iside = face information such as which points are on a face and which 
%                      which elements they belong to and, if a boundary, what
%                      is the boundary condition.  
%              jeside = for each element and each edge gives the FACE
%              number
%----------------------------------------------------------------------%
function [iside,jeside] = create_side(intma,bsido,npoin,nelem,nboun,nface,ngl)

%global arrays
iside = zeros(nface,4);
jeside= zeros(nelem,4);

%local arrays
lwher = zeros(npoin,1);
lhowm = zeros(npoin,1);
icone = zeros(5*npoin,1);
inode = zeros(4,1);
jnode = zeros(4,1);

%Fix lnode
inode(1)=1;
inode(2)=ngl;
inode(3)=ngl;
inode(4)=1;
jnode(1)=1;
jnode(2)=1;
jnode(3)=ngl;
jnode(4)=ngl;

%count how many elements own each node
for in=1:4
    for ie=1:nelem
        ip=intma(ie,inode(in),jnode(in));
        lhowm(ip)=lhowm(ip) + 1;
    end %ie
end %in

%track elements owning each node
lwher(1)=0;
for ip=2:npoin
   lwher(ip)=lwher(ip-1) + lhowm(ip-1);
end %ip

%another tracker array
lhowm = zeros(npoin,1);
for in=1:4
    for ie=1:nelem
        ip=intma(ie,inode(in),jnode(in));
        lhowm(ip)=lhowm(ip) + 1;
        jloca=lwher(ip) + lhowm(ip);
        icone(jloca)=ie;
    end %ie
end %in

%LOOP OVER THE NODES
iloca=0;
for ip=1:npoin
    iloc1=iloca;
    iele=lhowm(ip);
    
    if (iele ~= 0 ) 
        iwher=lwher(ip);

       %LOOP OVER THOSE ELEMENTS SURROUNDING NODE IP

       ip1=ip;
       for iel=1:iele
           ie=icone(iwher+iel);

           %find out position of ip in intma
           for in=1:4
               in1=in;
               ipt=intma(ie,inode(in),jnode(in));
               if (ipt == ip) 
                  break
               end
           end %in  
           
           %Check Edge of Element IE which claims IP
           j=0;
           for jnod=1:2:3
               iold=0;
               j=j+1;
               in2=in + jnod;
               if (in2 > 4) 
                  in2=in2-4;
               end
               ip2=intma(ie,inode(in2),jnode(in2));
               if (ip2 >= ip1) 

                  %check whether side is old or new
                  if (iloca ~= iloc1) 
                     for is=iloc1+1:iloca
                         iside(is,2);
                         jloca=is;
                         if (iside(is,2) == ip2) 
                            iold=1;
                            break;
                         end
                     end %is   
                  end %iloca   
                  
                  if (iold == 0)
                     %NEW SIDE
                     iloca=iloca + 1;
                     iside(iloca,1)=ip1;
                     iside(iloca,2)=ip2;
                     iside(iloca,2+j)=ie;
                  elseif (iold == 1)   
                     %OLD SIDE
                     iside(jloca,2+j)=ie;
                  end %iold
               end %ip2   
           end %jnod
       end %iel

       %Perform some Shifting to order the nodes of a side in CCW direction
       for is=iloc1+1:iloca
           if (iside(is,3) == 0) 
              iside(is,3)=iside(is,4);
              iside(is,4)=0;
              iside(is,1)=iside(is,2);
              iside(is,2)=ip1;
           end %iside   
       end %is
    end %if iele    
end %ip

if (iloca ~= nface) 
    disp( 'Error in SIDE. iloca nface = ');
    iloca
    nface
    pause
end 

%RESET THE BOUNDARY MARKERS
for is=1:nface
    if (iside(is,4) == 0) 
       il=iside(is,1);
       ir=iside(is,2);
       ie=iside(is,3);
       for ib=1:nboun
           ibe=bsido(ib,3);
           ibc=bsido(ib,4);
           if (ibe == ie) 
              ilb=bsido(ib,1);
              irb=bsido(ib,2);
              if  (ilb == il && irb == ir)
                  iside(is,4)=-ibc;
                  break
              end %ilb
           end %ibe
       end %ib
    end %iside
end %is

%FORM ELEMENT/SIDE CONNECTIVITY ARRAY
for is=1:nface
    iel=iside(is,3);
    ier=iside(is,4);
    is1=iside(is,1);
    is2=iside(is,2);

    %LEFT SIDE
    for in=1:4
       i1=intma(iel,inode(in),jnode(in));
       in1=in + 1;
       if (in1 > 4) 
          in1=1;
       end %in1
       i2=intma(iel,inode(in1),jnode(in1));
       if ((is1 == i1) && (is2 == i2)) 
          jeside(iel,in)=is;
       end %is1
    end %in   

    %RIGHT SIDE
    if (ier > 0) 
       for in=1:4
           i1=intma(ier,inode(in),jnode(in));
           in1=in + 1;
           if (in1 > 4) 
              in1=1;
           end %in1   
           i2=intma(ier,inode(in1),jnode(in1));
           if ((is1 == i2) && (is2 == i1)) 
              jeside(ier,in)=is;
           end %is1   
       end %in   
    end %ier
end %is

