%----------------------------------------------------------------------!
%This subroutine builds Periodic BCs along the 4 edges of a rectangular domain
%Written by Francis X. Giraldo on 2/2007
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%
% INPUT LIST: iside = face information
%             face = more face information
%             coord = gridpoint coordinates
%             nface = number of faces
%             nboun = number of boundaries
%
% OUTPUT LIST: face = augments the FACE data structure to include
%                     periodicity
%

%----------------------------------------------------------------------!
function face = create_face_periodicity(iside,face,coord,nface,nboun)

%Constant
tol=1e-6;

%Local arrays
ileft=zeros(nboun/4,1);
iright=zeros(nboun/4,1);
itop=zeros(nboun/4,1);
ibot=zeros(nboun/4,1);

%initialize
nleft=0; nright=0; ntop=0; nbot=0;

%Find Extrema of Domain
xmax=max(coord(:,1)); xmin=min(coord(:,1));
ymax=max(coord(:,2)); ymin=min(coord(:,2));

%loop thru sides and extract Left, Right, Bot, and Top
for is=1:nface

    %Check for Periodicity Edges
    ier=face(is,4);
    if (ier == -6) 

        i1=iside(is,1); i2=iside(is,2);
        xm=0.5*( coord(i1,1) + coord(i2,1) );
        ym=0.5*( coord(i1,2) + coord(i2,2) );

        %check Grid Point
        if ( abs(xm - xmin) < tol ) %left boundary
           nleft=nleft + 1;
           ileft(nleft)=is;
        elseif ( abs(xm - xmax) < tol ) %right boundary
           nright=nright + 1;
           iright(nright)=is;
        elseif ( abs(ym - ymin) < tol ) %bottom boundary
           nbot=nbot + 1;
           ibot(nbot)=is;
        elseif ( abs(ym - ymax) < tol ) %top boundary
           ntop=ntop + 1;
           itop(ntop)=is;
        else
           disp('No match in PERIODIC_BCS for is ier = ');
           is
           ier
           pause
        end %if
    end %ier
end %is

%Loop through Periodic BCs

%First: Do Left and Right
for i=1:nleft
    isl=ileft(i);
    i1=iside(isl,1);
    yl1=coord(i1,2);

    %Search for Corresponding Right Edge
    for j=1:nright
        isr=iright(j);
        i2=iside(isr,2);
        yr2=coord(i2,2);
        if ( abs(yl1-yr2) < tol ) %they match
           face(isl,2)=face(isr,1);
           face(isl,4)=face(isr,3);
           face(isr,3)=-6; %means skip it due to Periodicity
           iside(isl,4)=iside(isr,3);
           iside(isr,3)=-6;
           break;
        end %if
    end %j   
end %i

%Second: Do Top and Bottom
for i=1:nbot
    isl=ibot(i);
    i1=iside(isl,1);
    xl1=coord(i1,1);

    %Search for Corresponding Top Edge
    for j=1:ntop
        isr=itop(j);
        i2=iside(isr,2);
        xr2=coord(i2,1);
        if ( abs(xl1-xr2) < tol ) %they match
           face(isl,2)=face(isr,1);
           face(isl,4)=face(isr,3);
           face(isr,3)=-6; %means skip it due to Periodicity
           iside(isl,4)=iside(isr,3);
           iside(isr,3)=-6;
           break;
        end %if
     end %j   
end %i


