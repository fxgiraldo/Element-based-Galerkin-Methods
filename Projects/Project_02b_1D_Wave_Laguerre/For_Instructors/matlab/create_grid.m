%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 5/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [coord,coord_cg,intma,jac,wnq,npoin] = create_grid(ngl,nelem,nelem_LGL,xgl,xgr,wgl,wgr,icase,xgr_scale,basis_method_LGL_Only,xmin,xmax)

intma=zeros(ngl(nelem),nelem);
jac=zeros(nelem,1);
wnq=zeros(ngl(nelem),nelem);
coord=zeros(ngl(nelem),nelem);

%Set some constants
% if (icase == 1)
%     xmin=0;
%     xmax=50;
% end

% if basis_method_LGL_Only == 1
%     dx=(xmax-xmin)/(nelem);
% else
%     dx=(xmax-xmin)/(nelem-1);
% end
dx=(xmax-xmin)/(nelem_LGL);

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

      