%---------------------------------------------------------------------%
%This function computes the LGL grid and elements in 2D.
%Written by F.X. Giraldo on 4/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%
% INPUT LIST: nelx and nely are the number of elements in x and y
%             nop is the polynomial order
%             xgl are the interpolation points on the element.
% OUTPUT LIST:
%             coord are the coordinates: x=coord(:,1) and y=coord(:,2)
%             intma is the connectivity list that points to the global
%             gridpoint number
%             bsido is the boundary data (used by ISIDE and FACE)
%             iperiodic points to another point if periodicity is
%             applicable
%             npoin = number of global points
%             nelem = number of elements
%             nboun = number of boundary edges
%             nface = number of faces/edges in the grid
%---------------------------------------------------------------------%
function [coord,intma,bsido,iperiodic,npoin,nelem,nboun,nface] = create_grid_2d(nelx,nely,nop,xgl,plot_grid)

%Define Grid Dimensions
ngl=nop+1;
npoin=(nop*nelx + 1)*(nop*nely + 1);
nelem=nelx*nely;
nboun=2*nelx + 2*nely;
nface=2*nelem + nelx + nely;

%Initialize Global Arrays
coord=zeros(npoin,2);
intma=zeros(nelem,ngl,ngl);
bsido=zeros(nboun,4);
iperiodic=zeros(npoin,1);

%Initialize Local Arrays
node=zeros(npoin,npoin);

%Set some constants
dim=1.0;
xmin=-dim;
xmax=+dim;
ymin=-dim;
ymax=+dim;
dx=(xmax-xmin)/nelx;
dy=(ymax-ymin)/nely;
nop=ngl-1;
nx=nelx*nop + 1;
ny=nely*nop + 1;

%GENERATE COORD
ip=0;
jj=0;
for k=1:nely
   y0=ymin + real(k-1)*dy;

   if (k == 1) 
      l1=1;
   else
      l1=2;
   end

   for l=l1:ngl
      y=( xgl(l)+1 )*dy/2 + y0;

      jj=jj+1;
      ii=0;

      for i=1:nelx
         x0=xmin + real(i-1)*dx;

         if (i == 1)
            j1=1;
         else
            j1=2;
         end
   
         for j=j1:ngl
            ii=ii + 1;
            ip=ip + 1;
            x=( xgl(j)+1 )*dx/2 + x0;
            coord(ip,1)=x;
            coord(ip,2)=y;
            node(ii,jj)=ip;
         end %j
   
      end %i
   end %l   
end %k

%GENERATE INTMA
ie=0;
for k=1:nely
   for i=1:nelx
      ie=ie+1;
      for l=1:ngl
         jj=(ngl-1)*(k-1) + l;
         for j=1:ngl
            ii=(ngl-1)*(i-1) + j;
            ip=node(ii,jj);
            intma(ie,j,l)=ip;
         end %j
      end %l
   end %i   
end %k

%Generate BSIDO
ib=0;
for i=1:nelx
    ie=i;
    ib=ib+1;
    i1=(i-1)*(ngl-1) + 1;
    i2=(i-1)*(ngl-1) + ngl;
    ip1=node(i1,1);
    ip2=node(i2,1);
    bsido(ib,1)=ip1;
    bsido(ib,2)=ip2;
    bsido(ib,3)=ie;
    bsido(ib,4)=6;
end

%Right Boundary
for i=1:nely
    ie=(nelx)*(i);
    ib=ib+1;
    i1=(i-1)*(ngl-1) + 1;
    i2=(i-1)*(ngl-1) + ngl;
    ip1=node(nx,i1);
    ip2=node(nx,i2);
    bsido(ib,1)=ip1;
    bsido(ib,2)=ip2;
    bsido(ib,3)=ie;
    bsido(ib,4)=6;
end 

%Top Boundary
for i=nelx:-1:1 
    ie=nelem - (nelx - i);
    ib=ib+1;
    i1=(i-1)*(ngl-1) + ngl;
    i2=(i-1)*(ngl-1) + 1;
    ip1=node(i1,ny);
    ip2=node(i2,ny);
    bsido(ib,1)=ip1;
    bsido(ib,2)=ip2;
    bsido(ib,3)=ie;
    bsido(ib,4)=6;
end 

%Left Boundary
for i=nely:-1:1
    ie=(nelx)*(i-1) + 1;
    ib=ib+1;
    i1=(i-1)*(ngl-1) + ngl;
    i2=(i-1)*(ngl-1) + 1;
    ip1=node(1,i1);
    ip2=node(1,i2);
    bsido(ib,1)=ip1;
    bsido(ib,2)=ip2;
    bsido(ib,3)=ie;
    bsido(ib,4)=6;
end

%Periodicity
for i=1:npoin
   iperiodic(i)=i;
end

%X-Periodicity
for i=1:ny
   i1=node(1,i);
   i2=node(nx,i);
   iperiodic(i2)=i1;
end

%Y-Periodicity
for i=1:nx
   i1=node(i,1);
   i2=node(i,ny);
   iperiodic(i2)=iperiodic(i1);
end

%Plot Grid
if (plot_grid == 1)
    x=zeros(5,1);
    y=zeros(5,1);
    figure;
    hold on;
    for e=1:nelem
        for j=1:ngl-1
            for i=1:ngl-1
                i1=intma(e,i,j);
                i2=intma(e,i+1,j);
                i3=intma(e,i+1,j+1);
                i4=intma(e,i,j+1);
                x(1)=coord(i1,1); y(1)=coord(i1,2);
                x(2)=coord(i2,1); y(2)=coord(i2,2);
                x(3)=coord(i3,1); y(3)=coord(i3,2);
                x(4)=coord(i4,1); y(4)=coord(i4,2);
                x(5)=coord(i1,1); y(5)=coord(i1,2);
                plot_handle=plot(x,y,'-r');
                set(plot_handle,'LineWidth',1.5);
            end
        end
        i1=intma(e,1,1);
        i2=intma(e,ngl,1);
        i3=intma(e,ngl,ngl);
        i4=intma(e,1,ngl);
        x(1)=coord(i1,1); y(1)=coord(i1,2);
        x(2)=coord(i2,1); y(2)=coord(i2,2);
        x(3)=coord(i3,1); y(3)=coord(i3,2);
        x(4)=coord(i4,1); y(4)=coord(i4,2);
        x(5)=coord(i1,1); y(5)=coord(i1,2);
        plot_handle=plot(x,y,'-b');
        set(plot_handle,'LineWidth',2);
    end
    title_text=['Original Grid Plot For: Ne = ' num2str(nelem) ', N = ' num2str(nop)];
    title([title_text],'FontSize',18);      

    xlabel('X','FontSize',18);
    ylabel('Y','FontSize',18);
    axis image
end 


      


      
