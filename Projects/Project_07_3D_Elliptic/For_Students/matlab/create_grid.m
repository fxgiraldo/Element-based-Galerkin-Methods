%---------------------------------------------------------------------%
%This function computes the LGL grid and elements in 2D.
%Written by F.X. Giraldo on 4/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [coord,intma,iboun] = create_grid(npoin,nelem,nboun,nelx,nely,nelz,ngl,xgl,plot_grid,lwarp_grid)

%Initialize Global Arrays
coord=zeros(3,npoin);
intma=zeros(ngl,ngl,ngl,nelem);
iboun=zeros(nboun,1);

%Set some constants
xmin=-1;
xmax=+1;
ymin=-1;
ymax=+1;
zmin=-1;
zmax=+1;
dx=(xmax-xmin)/nelx;
dy=(ymax-ymin)/nely;
dz=(zmax-zmin)/nelz;
nop=ngl-1;
nx=nelx*nop + 1;
ny=nely*nop + 1;
nz=nelz*nop + 1;

%Initialize Local Arrays
node=zeros(nx,ny,nz);

%GENERATE COORD
I=0;
kk=0;
for m=1:nelz
    z0=zmin + real(m-1)*dz;
    if (m == 1)
        n1=1;
    else
        n1=2;
    end
    
    for n=n1:ngl
        z=( xgl(n)+1 )*dz/2 + z0;
        kk=kk+1;
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
                        I=I + 1;
                        x=( xgl(j)+1 )*dx/2 + x0;
                        coord(1,I)=x;
                        coord(2,I)=y;
                        coord(3,I)=z;
                        node(ii,jj,kk)=I;
                    end %j
                    
                end %i
            end %l
        end %k
    end %n
end %m

%GENERATE INTMA
e=0;
for m=1:nelz
    for k=1:nely
        for i=1:nelx
            e=e+1;
            for n=1:ngl
                kk=(ngl-1)*(m-1) + n;
                for l=1:ngl
                    jj=(ngl-1)*(k-1) + l;
                    for j=1:ngl
                        ii=(ngl-1)*(i-1) + j;
                        I=node(ii,jj,kk);
                        intma(j,l,n,e)=I;
                    end %j
                end %n
            end %l
        end %i
    end %k
end %m

%GENERATE Boundary Conditions
B=0;

%Front Face (y=-1)
j=1;
for k=1:nz
    for i=1:nx
        I=node(i,j,k);
        B=B + 1;
        iboun(B)=I;
    end
end

%Back Face (y=+1)
j=ny;
for k=1:nz
    for i=1:nx
        I=node(i,j,k);
        B=B + 1;
        iboun(B)=I;
    end
end


%Left Face (x=-1)
i=1;
for k=1:nz
    for j=2:ny-1
        I=node(i,j,k);
        B=B + 1;
        iboun(B)=I;
    end
end

%Right Face (x=+1)
i=nx;
for k=1:nz
    for j=2:ny-1
        I=node(i,j,k);
        B=B + 1;
        iboun(B)=I;
    end
end

%Bottom Face (z=-1)
k=1;
for j=2:ny-1
    for i=2:nx-1
        I=node(i,j,k);
        B=B + 1;
        iboun(B)=I;
    end
end

%Top Face (z=-1)
k=nz;
for j=2:ny-1
    for i=2:nx-1
        I=node(i,j,k);
        B=B + 1;
        iboun(B)=I;
    end
end

if (B ~= nboun) 
    disp('Error in Create_Grid');
    return
end

%Warp Grid
if (lwarp_grid == 1)
    coord = warp_grid(coord,npoin);
end

%Plot Grid
if (plot_grid == 1)
    x=zeros(5,1);
    y=zeros(5,1);
    z=zeros(5,1);
    figure;
    hold on;
    for e=1:nelem
        
        %---------------------Build Lobatto Points------------------------%
        %Build x-y plane grid for each z value
        for k=1:ngl
            for j=1:ngl-1
                for i=1:ngl-1
                    i1=intma(i,j,k,e);
                    i2=intma(i+1,j,k,e);
                    i3=intma(i+1,j+1,k,e);
                    i4=intma(i,j+1,k,e);
                    x(1)=coord(1,i1); y(1)=coord(2,i1); z(1)=coord(3,i1);
                    x(2)=coord(1,i2); y(2)=coord(2,i2); z(2)=coord(3,i2);
                    x(3)=coord(1,i3); y(3)=coord(2,i3); z(3)=coord(3,i3);
                    x(4)=coord(1,i4); y(4)=coord(2,i4); z(4)=coord(3,i4);
                    x(5)=coord(1,i1); y(5)=coord(2,i1); z(5)=coord(3,i1);
                    plot_handle=plot3(x,y,z,':r');
                    set(plot_handle,'LineWidth',1);
                end
            end
        end

        %Build x-z plane grid for each y value
        for j=1:ngl
            for k=1:ngl-1
                for i=1:ngl-1
                    i1=intma(i,j,k,e);
                    i2=intma(i+1,j,k,e);
                    i3=intma(i+1,j,k+1,e);
                    i4=intma(i,j,k+1,e);
                    x(1)=coord(1,i1); y(1)=coord(2,i1); z(1)=coord(3,i1);
                    x(2)=coord(1,i2); y(2)=coord(2,i2); z(2)=coord(3,i2);
                    x(3)=coord(1,i3); y(3)=coord(2,i3); z(3)=coord(3,i3);
                    x(4)=coord(1,i4); y(4)=coord(2,i4); z(4)=coord(3,i4);
                    x(5)=coord(1,i1); y(5)=coord(2,i1); z(5)=coord(3,i1);
                    plot_handle=plot3(x,y,z,':r');
                    set(plot_handle,'LineWidth',1);
                end
            end
        end

        %Build y-z plane grid for each x value
        for i=1:ngl
            for k=1:ngl-1
                for j=1:ngl-1
                    i1=intma(i,j,k,e);
                    i2=intma(i,j+1,k,e);
                    i3=intma(i,j+1,k+1,e);
                    i4=intma(i,j,k+1,e);
                    x(1)=coord(1,i1); y(1)=coord(2,i1); z(1)=coord(3,i1);
                    x(2)=coord(1,i2); y(2)=coord(2,i2); z(2)=coord(3,i2);
                    x(3)=coord(1,i3); y(3)=coord(2,i3); z(3)=coord(3,i3);
                    x(4)=coord(1,i4); y(4)=coord(2,i4); z(4)=coord(3,i4);
                    x(5)=coord(1,i1); y(5)=coord(2,i1); z(5)=coord(3,i1);
                    plot_handle=plot3(x,y,z,':r');
                    set(plot_handle,'LineWidth',1);
                end
            end
        end

        %---------------------Build Element Boundary----------------------%
        %Build x-y plane grid for k=1 value
        i1=intma(1,1,1,e);
        i2=intma(ngl,1,1,e);
        i3=intma(ngl,ngl,1,e);
        i4=intma(1,ngl,1,e);
        x(1)=coord(1,i1); y(1)=coord(2,i1); z(1)=coord(3,i1);
        x(2)=coord(1,i2); y(2)=coord(2,i2); z(2)=coord(3,i2);
        x(3)=coord(1,i3); y(3)=coord(2,i3); z(3)=coord(3,i3);
        x(4)=coord(1,i4); y(4)=coord(2,i4); z(4)=coord(3,i4);
        x(5)=coord(1,i1); y(5)=coord(2,i1); z(5)=coord(3,i1);
        plot_handle=plot3(x,y,z,'-b');
        set(plot_handle,'LineWidth',1);

        %Build x-y plane grid for k=ngl value
        i1=intma(1,1,ngl,e);
        i2=intma(ngl,1,ngl,e);
        i3=intma(ngl,ngl,ngl,e);
        i4=intma(1,ngl,ngl,e);
        x(1)=coord(1,i1); y(1)=coord(2,i1); z(1)=coord(3,i1);
        x(2)=coord(1,i2); y(2)=coord(2,i2); z(2)=coord(3,i2);
        x(3)=coord(1,i3); y(3)=coord(2,i3); z(3)=coord(3,i3);
        x(4)=coord(1,i4); y(4)=coord(2,i4); z(4)=coord(3,i4);
        x(5)=coord(1,i1); y(5)=coord(2,i1); z(5)=coord(3,i1);
        plot_handle=plot3(x,y,z,'-b');
        set(plot_handle,'LineWidth',1);

        %Build x-z plane grid for j=1
        i1=intma(1,1,1,e);
        i2=intma(ngl,1,1,e);
        i3=intma(ngl,1,ngl,e);
        i4=intma(1,1,ngl,e);
        x(1)=coord(1,i1); y(1)=coord(2,i1); z(1)=coord(3,i1);
        x(2)=coord(1,i2); y(2)=coord(2,i2); z(2)=coord(3,i2);
        x(3)=coord(1,i3); y(3)=coord(2,i3); z(3)=coord(3,i3);
        x(4)=coord(1,i4); y(4)=coord(2,i4); z(4)=coord(3,i4);
        x(5)=coord(1,i1); y(5)=coord(2,i1); z(5)=coord(3,i1);
        plot_handle=plot3(x,y,z,'-b');
        set(plot_handle,'LineWidth',1);

        %Build x-z plane grid for j=ngl
        i1=intma(1,ngl,1,e);
        i2=intma(ngl,ngl,1,e);
        i3=intma(ngl,ngl,ngl,e);
        i4=intma(1,ngl,ngl,e);
        x(1)=coord(1,i1); y(1)=coord(2,i1); z(1)=coord(3,i1);
        x(2)=coord(1,i2); y(2)=coord(2,i2); z(2)=coord(3,i2);
        x(3)=coord(1,i3); y(3)=coord(2,i3); z(3)=coord(3,i3);
        x(4)=coord(1,i4); y(4)=coord(2,i4); z(4)=coord(3,i4);
        x(5)=coord(1,i1); y(5)=coord(2,i1); z(5)=coord(3,i1);
        plot_handle=plot3(x,y,z,'-b');
        set(plot_handle,'LineWidth',1);

        %Build y-z plane grid for each i=1 value
        i1=intma(1,1,1,e);
        i2=intma(1,ngl,1,e);
        i3=intma(1,ngl,ngl,e);
        i4=intma(1,1,ngl,e);
        x(1)=coord(1,i1); y(1)=coord(2,i1); z(1)=coord(3,i1);
        x(2)=coord(1,i2); y(2)=coord(2,i2); z(2)=coord(3,i2);
        x(3)=coord(1,i3); y(3)=coord(2,i3); z(3)=coord(3,i3);
        x(4)=coord(1,i4); y(4)=coord(2,i4); z(4)=coord(3,i4);
        x(5)=coord(1,i1); y(5)=coord(2,i1); z(5)=coord(3,i1);
        plot_handle=plot3(x,y,z,'-b');
        set(plot_handle,'LineWidth',1);

        %Build y-z plane grid for each i=ngl value
        i1=intma(ngl,1,1,e);
        i2=intma(ngl,ngl,1,e);
        i3=intma(ngl,ngl,ngl,e);
        i4=intma(ngl,1,ngl,e);
        x(1)=coord(1,i1); y(1)=coord(2,i1); z(1)=coord(3,i1);
        x(2)=coord(1,i2); y(2)=coord(2,i2); z(2)=coord(3,i2);
        x(3)=coord(1,i3); y(3)=coord(2,i3); z(3)=coord(3,i3);
        x(4)=coord(1,i4); y(4)=coord(2,i4); z(4)=coord(3,i4);
        x(5)=coord(1,i1); y(5)=coord(2,i1); z(5)=coord(3,i1);
        plot_handle=plot3(x,y,z,'-b');
        set(plot_handle,'LineWidth',1);
    end
    title_text=['Grid Plot For: Ne = ' num2str(nelem) ', N = ' num2str(nop)];
    title([title_text],'FontSize',18);      

    xlabel('X','FontSize',18);
    ylabel('Y','FontSize',18);
    zlabel('Z','FontSize',18);
    axis image
end 
