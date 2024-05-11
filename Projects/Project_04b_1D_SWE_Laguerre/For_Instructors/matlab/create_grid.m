%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 5/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [coord,coord_cg,intma,jac,wnq,npoin,xmin,xmax] = create_grid(ngl,ngl_LGR,nelem,nelem_LGL,xgl,xgr,wgl,wgr,icase,LGR_basis,LGR_scale,LGR_scale_factor,LGR_artificial_layer,lsponge,delta_nl)

intma=zeros(ngl(nelem),nelem);
jac=zeros(nelem,1);
wnq=zeros(ngl(nelem),nelem);
coord=zeros(ngl(nelem),nelem);

%Set some constants
basis_method_LGL_Only=0;
if strcmp(LGR_basis,'LGL')
    basis_method_LGL_Only=1;
end

xmin=-1;
xmax=+1;
if (icase == 0)
    xmin=-1;
    xmax=+1;
elseif (icase == 1)
    if delta_nl==0
        xmin=-2;
        xmax=+1;
    elseif delta_nl==1
        xmin=-2;
        xmax=+1;
    end
elseif (icase == 2)
    xmin=-12;
    xmax=+3;
elseif (icase == 3)
    xmin=0;
    xmax=1;
elseif (icase == 4)
    xmin=-3;
    xmax=+3;
end
%Scale LGR element to be closer to LGL elements
if LGR_artificial_layer==1
    dx_lgl=(xmax-xmin)/nelem_LGL;
else
    dx_lgl=(xmax-xmin)/nelem;
end
dx_lgr=xgr(ngl_LGR)-xgr(1);
xgr_scale=1.0;
if (LGR_scale==1)
    xgr_scale=dx_lgl/dx_lgr*LGR_scale_factor;
end

if (lsponge==0)
    dx=(xmax-xmin)/nelem;
elseif (lsponge==1)
    dx=dx_lgl;
end

%Generate LGL Grid Points
I=1;
coord(1,1)=xmin;
for i=1:nelem_LGL
   x0=xmin + (i-1)*dx;
   coord(1,i)=x0;
   intma(1,i)=I;
   jac(i)=dx/2;
   for j=2:ngl(i)
      I=I+1;
      coord(j,i)=x0 + ( xgl(j)+1 )*dx/2;
      intma(j,i)=I;
   end %j
end %i

%Generate LGR Grid Points
dx=dx_lgr*xgr_scale;
for i=nelem_LGL+1:nelem
    jac(i)=(1-basis_method_LGL_Only)*xgr_scale + basis_method_LGL_Only*dx/2;
    x0=coord(ngl(i-1),i-1);
    coord(1,i)=x0;
    intma(1,i)=I;
    for j=2:ngl(i)
        I=I + 1;
        coord(j,i)=x0 + (1-basis_method_LGL_Only)*xgr_scale*xgr(j) + basis_method_LGL_Only*( xgr(j)+1 )*dx/2;
        intma(j,i)=I;
    end %j
end %i

npoin=I;

%Generate CG coord
coord_cg=zeros(npoin,1);
for e=1:nelem
    for i=1:ngl(e)
        I=intma(i,e);
        coord_cg(I)=coord(i,e);
    end
end

%Form Unified Quadrature Weights
for i=1:nelem_LGL
    for j=1:ngl(i)
        wnq(j,i)=wgl(j);
    end
end
for i=nelem_LGL+1:nelem
    for j=1:ngl(i)
        wnq(j,i)=wgr(j);
    end
end