%---------------------------------------------------------------------%
%This function computes the LGL grid and elements in 2D.
%Written by F.X. Giraldo on 4/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [coord,intma,bsido,iperiodic] = create_grid_2d(npoin,nelem,nboun,nelx,nely,ngl,xgl,eq_set)

%Initialize Global Arrays
coord=zeros(2,npoin);
intma=zeros(ngl,ngl,nelem);
bsido=zeros(4,nboun);
iperiodic=zeros(npoin,1);

%Initialize Local Arrays
node=zeros(npoin,npoin);

boun=4; %no-flux
if strcmp(eq_set,'advection')
   boun=6; %periodic         
end

%Set some constants
xmin=-1;
xmax=+1;
ymin=-1;
ymax=+1;
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
            coord(1,ip)=x;
            coord(2,ip)=y;
            node(ii,jj)=ip;
         end %j
   
      end %i
   end %l   
end %k

%GENERATE INTMA
e=0;
for k=1:nely
   for i=1:nelx
      e=e+1;
      for l=1:ngl
         jj=(ngl-1)*(k-1) + l;
         for j=1:ngl
            ii=(ngl-1)*(i-1) + j;
            ip=node(ii,jj);
            intma(j,l,e)=ip;
         end %j
      end %l
   end %i   
end %k

%Generate BSIDO
ib=0;
for i=1:nelx
    e=i;
    ib=ib+1;
    i1=(i-1)*(ngl-1) + 1;
    i2=(i-1)*(ngl-1) + ngl;
    ip1=node(i1,1);
    ip2=node(i2,1);
    bsido(1,ib)=ip1;
    bsido(2,ib)=ip2;
    bsido(3,ib)=e;
    bsido(4,ib)=boun;
end

%Right Boundary
for i=1:nely
    e=(nelx)*(i);
    ib=ib+1;
    i1=(i-1)*(ngl-1) + 1;
    i2=(i-1)*(ngl-1) + ngl;
    ip1=node(nx,i1);
    ip2=node(nx,i2);
    bsido(1,ib)=ip1;
    bsido(2,ib)=ip2;
    bsido(3,ib)=e;
    bsido(4,ib)=boun;
end 

%Top Boundary
for i=nelx:-1:1 
    e=nelem - (nelx - i);
    ib=ib+1;
    i1=(i-1)*(ngl-1) + ngl;
    i2=(i-1)*(ngl-1) + 1;
    ip1=node(i1,ny);
    ip2=node(i2,ny);
    bsido(1,ib)=ip1;
    bsido(2,ib)=ip2;
    bsido(3,ib)=e;
    bsido(4,ib)=boun;
end 

%Left Boundary
for i=nely:-1:1
    e=(nelx)*(i-1) + 1;
    ib=ib+1;
    i1=(i-1)*(ngl-1) + ngl;
    i2=(i-1)*(ngl-1) + 1;
    ip1=node(1,i1);
    ip2=node(1,i2);
    bsido(1,ib)=ip1;
    bsido(2,ib)=ip2;
    bsido(3,ib)=e;
    bsido(4,ib)=boun;
end

%Periodicity
for i=1:npoin
   iperiodic(i)=i;
end

if strcmp(eq_set,'advection')
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
end    