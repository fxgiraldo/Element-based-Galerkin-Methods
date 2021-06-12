%---------------------------------------------------------------------%
%This function computes the LGL grid and elements in 2D.
%Written by F.X. Giraldo on 4/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [coord,intma,iboun,iperiodic] = create_grid(npoin,nelem,nboun,nelx,nely,ngl,xgl,plot_grid,rotate_grid)

%Initialize Global Arrays
coord=zeros(2,npoin);
intma=zeros(ngl,ngl,nelem);
iboun=zeros(nboun,1);
iperiodic=zeros(npoin,1);

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

%GENERATE Boundary Conditions
ib=0;
for i=1:nx
    ip=node(i,1);
    ib=ib + 1;
    iboun(ib)=ip;
end
for j=2:ny
    ip=node(nx,j);
    ib=ib + 1;
    iboun(ib)=ip;
end
for i=nx-1:-1:1
    ip=node(i,ny);
    ib=ib + 1;
    iboun(ib)=ip;
end
for j=ny-1:-1:2
    ip=node(1,j);
    ib=ib + 1;
    iboun(ib)=ip;
end

if (ib ~= nboun) 
    disp('Error in Create_Grid');
    return
end

%Periodicity
for i=1:npoin
   iperiodic(i)=i;
end

%X-Periodicity
for i=1:ny
   i1=node(1,i);
   i2=node(nx,i);
   %iperiodic(i2)=i1;
end

%Y-Periodicity
for i=1:nx
   i1=node(i,1);
   i2=node(i,ny);
   %iperiodic(i2)=i1;
end

%Rotate Grid
if (rotate_grid == 1)
    coord = rotate_grid_function(coord,npoin,ngl);
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
                i1=intma(i,j,e);
                i2=intma(i+1,j,e);
                i3=intma(i+1,j+1,e);
                i4=intma(i,j+1,e);
                x(1)=coord(1,i1); y(1)=coord(2,i1);
                x(2)=coord(1,i2); y(2)=coord(2,i2);
                x(3)=coord(1,i3); y(3)=coord(2,i3);
                x(4)=coord(1,i4); y(4)=coord(2,i4);
                x(5)=coord(1,i1); y(5)=coord(2,i1);
                plot_handle=plot(x,y,'-r');
                set(plot_handle,'LineWidth',1);
            end
        end
        i1=intma(1,1,e);
        i2=intma(ngl,1,e);
        i3=intma(ngl,ngl,e);
        i4=intma(1,ngl,e);
        x(1)=coord(1,i1); y(1)=coord(2,i1);
        x(2)=coord(1,i2); y(2)=coord(2,i2);
        x(3)=coord(1,i3); y(3)=coord(2,i3);
        x(4)=coord(1,i4); y(4)=coord(2,i4);
        x(5)=coord(1,i1); y(5)=coord(2,i1);
        plot_handle=plot(x,y,'-b');
        set(plot_handle,'LineWidth',2);
    end
    title_text=['Grid Plot For: Ne = ' num2str(nelem) ', N = ' num2str(nop)];
    title([title_text],'FontSize',18);      

    xlabel('X','FontSize',18);
    ylabel('Y','FontSize',18);
    axis image
end 
