%----------------------------------------------------------------------!
%This subroutine builds Periodic BCs along the 4 edges of a rectangular domain
%Written by Francis X. Giraldo on 2/2007
%           Naval Postgraduate School
%           Department of Applied Mathematics
%           Monterey, CA 93943-5216
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
xmax=max(coord(1,:)); xmin=min(coord(1,:));
ymax=max(coord(2,:)); ymin=min(coord(2,:));

%loop thru sides and extract Left, Right, Bot, and Top
for is=1:nface

    %Check for Periodicity Edges
    ier=face(4,is);
    if (ier == -6) 

        i1=iside(1,is); i2=iside(2,is);
        xm=0.5*( coord(1,i1) + coord(1,i2) );
        ym=0.5*( coord(2,i1) + coord(2,i2) );

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
    i1=iside(1,isl);
    yl1=coord(2,i1);

    %Search for Corresponding Right Edge
    for j=1:nright
        isr=iright(j);
        i2=iside(2,isr);
        yr2=coord(2,i2);
        if ( abs(yl1-yr2) < tol ) %they match
           face(2,isl)=face(1,isr);
           face(4,isl)=face(3,isr);
           face(3,isr)=-6; %means skip it due to Periodicity
           iside(4,isl)=iside(3,isr);
           iside(3,isr)=-6;
           break;
        end %if
    end %j   
end %i

%Second: Do Top and Bottom
for i=1:nbot
    isl=ibot(i);
    i1=iside(1,isl);
    xl1=coord(1,i1);

    %Search for Corresponding Top Edge
    for j=1:ntop
        isr=itop(j);
        i2=iside(2,isr);
        xr2=coord(1,i2);
        if ( abs(xl1-xr2) < tol ) %they match
           face(2,isl)=face(1,isr);
           face(4,isl)=face(3,isr);
           face(3,isr)=-6; %means skip it due to Periodicity
           iside(4,isl)=iside(3,isr);
           iside(3,isr)=-6;
           break;
        end %if
     end %j   
end %i


