%---------------------------------------------------------------------%
%This function computes the LGL grid and elements in 2D.
%Written by F.X. Giraldo on 4/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [coord,intma,bsido,iperiodic,bsidoa] = create_grid_2d_dg(npoin,nelem,nboun,nelx,nely,ngl,xgl,ae)

%Initialize Global Arrays - four times as long as without refinement
coord=zeros(npoin*ae,2);
intma=zeros(nelem*ae,ngl,ngl);
bsido=zeros(nboun,ae);
bsidoa=zeros(nboun/2,5);
iperiodic=zeros(npoin*ae,1);

%Initialize Local Arrays
node=zeros(npoin,npoin);
ibc=zeros(npoin);

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
    
    %---
    bsidoa(ib,1) = 3;
    bsidoa(ib,2) = 1;
    bsidoa(ib,3) = ie+nelx*(nely-1);
    bsidoa(ib,4) = ie;
    bsidoa(ib,5) = 6;
    %---
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
    
    %---
    bsidoa(ib,1) = 2;
    bsidoa(ib,2) = 4;
    bsidoa(ib,3) = ie;
    bsidoa(ib,4) = ie-(nely-1);
    bsidoa(ib,5) = 6;
    %---
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
   iperiodic(i2)=i1;
end



      


      
