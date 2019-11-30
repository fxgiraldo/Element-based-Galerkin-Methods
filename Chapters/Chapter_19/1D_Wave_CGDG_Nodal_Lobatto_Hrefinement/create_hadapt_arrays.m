%---------------------------------------------------------------------%
%This function creates the H-refinement arrays for all levels of
%refinement.
%Written by F.X. Giraldo on 3/10/2016
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [parent,children,active,ref_level,elem_lev_pointer,coord,intma,nelem,nelem0] = create_hadapt_arrays(intma,coord,nelem,ngl,xgl,hadapt_lev)

%Initialize values
elem_lev_pointer=zeros(2,hadapt_lev+1);
elem_lev_pointer(1,1)=1;
elem_lev_pointer(2,1)=nelem;
ii=1;
for i=1:hadapt_lev
    ii=ii + 2^(i);
    elem_lev_pointer(1,i+1)=elem_lev_pointer(2,i)+1;
    elem_lev_pointer(2,i+1)=elem_lev_pointer(2,i) + 2^(i)*nelem;
end
nelem_lev=nelem*ii;

parent=zeros(nelem_lev,1);
children=zeros(2,nelem_lev);
active=zeros(nelem_lev,1);
coord_lev=zeros(ngl,nelem_lev);
intma_lev=zeros(ngl,nelem_lev);
ref_level=zeros(nelem_lev,1);

%Initialize Grid
for e=1:nelem
    coord_lev(:,e)=coord(:,e);
    intma_lev(:,e)=intma(:,e);
    active(e)=1;
    ref_level(e)=0;
end

%Generate Refined grids
ip=nelem*ngl;
for iloop=1:hadapt_lev
    ii=1;
    for i=1:iloop-1
        ii=ii + 2^(i);
    end
    jj=1;
    for i=1:iloop
        jj=jj + 2^(i);
    end
    iloop;
    nelem_jj=nelem*jj;
    nelem_new=nelem*ii;
    
    for e=elem_lev_pointer(1,iloop):elem_lev_pointer(2,iloop)
        %Form New Coords
        x0=coord_lev(1,e);
        xN=coord_lev(ngl,e);
        xm=0.5*(x0+xN);
        dx=0.5*(xN-x0);
        
        %Left Element
        nelem_left=nelem_new+1;
        ref_level(nelem_left)=iloop;
        for j=1:ngl
            coord_lev(j,nelem_left)=( xgl(j)+1 )*dx/2 + x0;
            ip=ip+1;
            intma_lev(j,nelem_left)=ip;
        end
        %Right Element
        nelem_right=nelem_new+2;
        ref_level(nelem_right)=iloop;
        for j=1:ngl
            coord_lev(j,nelem_right)=( xgl(j)+1 )*dx/2 + xm;
            ip=ip+1;
            intma_lev(j,nelem_right)=ip;
        end
        
        %Create Tree Arrays
        parent(nelem_left)=e; parent(nelem_right)=e;
        children(1,e)=nelem_left; children(2,e)=nelem_right;
        nelem_new=nelem_new + 2;
    end %e
    %check if size is still correct
    nelem_new;
    nelem_jj;
    if (nelem_new ~= nelem_jj)
        disp(['Error in CREATE_HDADAPT_ARRAYS: nelem_new = ',num2str(nelem_new),' should be = ',num2str(nelem_jj)]);
        pause;
    end
end %iloop

%Rewrite Grid arrays
nelem0=nelem;
nelem=nelem_lev;
coord=coord_lev;
intma=intma_lev;