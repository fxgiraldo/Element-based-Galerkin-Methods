%----------------------------------------------------------------------%
%This subroutine creates the array ISIDE which stores all of
%the information concerning the sides of all the elements.
%Written by F.X. Giraldo on 5/2021
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%----------------------------------------------------------------------%
function [iside,jeside] = create_side(intma,bsido,npoin,nelem,nboun,nside,ngl)

%global arrays
iside = zeros(4,nside);
jeside= zeros(4,nelem);

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
    for e=1:nelem
        I=intma(inode(in),jnode(in),e);
        lhowm(I)=lhowm(I) + 1;
    end %ie
end %in

%track elements owning each node
lwher(1)=0;
for I=2:npoin
   lwher(I)=lwher(I-1) + lhowm(I-1);
end %ip

%another tracker array
lhowm = zeros(npoin,1);
for in=1:4
    for e=1:nelem
        I=intma(inode(in),jnode(in),e);
        lhowm(I)=lhowm(I) + 1;
        jloca=lwher(I) + lhowm(I);
        icone(jloca)=e;
    end %ie
end %in

%LOOP OVER THE NODES
iloca=0;
for I=1:npoin
    iloc1=iloca;
    iele=lhowm(I);
    
    if (iele ~= 0 ) 
        iwher=lwher(I);

       %LOOP OVER THOSE ELEMENTS SURROUNDING NODE I

       ip1=I;
       for el=1:iele
           e=icone(iwher+el);

           %find out position of ip in intma
           for in=1:4
               ipt=intma(inode(in),jnode(in),e);
               if (ipt == I) 
                  break
               end
           end %in  
           
           %Check Edge of Element E which claims I
           j=0;
           for jnod=1:2:3
               iold=0;
               j=j+1;
               in2=in + jnod;
               if (in2 > 4) 
                  in2=in2-4;
               end
               ip2=intma(inode(in2),jnode(in2),e);
               if (ip2 >= ip1) 

                  %check whether side is old or new
                  if (iloca ~= iloc1) 
                     for is=iloc1+1:iloca
                         jloca=is;
                         if (iside(2,is) == ip2) 
                            iold=1;
                            break;
                         end
                     end %is   
                  end %iloca   
                  
                  if (iold == 0)
                     %NEW SIDE
                     iloca=iloca + 1;
                     iside(1,iloca)=ip1;
                     iside(2,iloca)=ip2;
                     iside(2+j,iloca)=e;
                  elseif (iold == 1)   
                     %OLD SIDE
                     iside(2+j,jloca)=e;
                  end %iold
               end %ip2   
           end %jnod
       end %iel

       %Perform some Shifting to order the nodes of a side in CCW direction
       for is=iloc1+1:iloca
           if (iside(3,is) == 0) 
              iside(3,is)=iside(4,is);
              iside(4,is)=0;
              iside(1,is)=iside(2,is);
              iside(2,is)=ip1;
           end %iside   
       end %is
    end %if iele    
end %ip

if (iloca ~= nside) 
    disp( 'Error in SIDE. iloca nside = ');
    iloca
    nside
    pause
end 

%RESET THE BOUNDARY MARKERS
for is=1:nside
    if (iside(4,is) == 0) 
       il=iside(1,is);
       ir=iside(2,is);
       e=iside(3,is);
       for ib=1:nboun
           ibe=bsido(3,ib);
           ibc=bsido(4,ib);
           if (ibe == e) 
              ilb=bsido(1,ib);
              irb=bsido(2,ib);
              if  (ilb == il && irb == ir)
                  iside(4,is)=-ibc;
                  break
              end %ilb
           end %ibe
       end %ib
    end %iside
end %is

%FORM ELEMENT/SIDE CONNECTIVITY ARRAY
for is=1:nside
    el=iside(3,is);
    er=iside(4,is);
    is1=iside(1,is);
    is2=iside(2,is);

    %LEFT SIDE
    for in=1:4
       i1=intma(inode(in),jnode(in),el);
       in1=in + 1;
       if (in1 > 4) 
          in1=1;
       end %in1
       i2=intma(inode(in1),jnode(in1),el);
       if ((is1 == i1) && (is2 == i2)) 
          jeside(in,el)=is;
       end %is1
    end %in   

    %RIGHT SIDE
    if (er > 0) 
       for in=1:4
           i1=intma(inode(in),jnode(in),er);
           in1=in + 1;
           if (in1 > 4) 
              in1=1;
           end %in1   
           i2=intma(inode(in1),jnode(in1),er);
           if ((is1 == i2) && (is2 == i1)) 
              jeside(in,er)=is;
           end %is1   
       end %in   
    end %ier
end %is

