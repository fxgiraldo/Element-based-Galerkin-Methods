%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 5/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [coord_cg,coord_dg,intma_cg,intma_dg,periodicity_cg,periodicity_dg] = create_grid(ngl,nelem,npoin_cg,npoin_dg,xgl)

%Set some constants
xmin=-1;
xmax=+1;
dx=(xmax-xmin)/nelem;
coord_cg=zeros(npoin_cg,1);
coord_dg=zeros(npoin_dg,1);

%Generate COORD and INTMA for CG
ip=1;
coord_cg(1)=xmin;
for e=1:nelem
   x0=xmin + (e-1)*dx;
   intma_cg(1,e)=ip;
   for i=2:ngl
      ip=ip + 1;
      coord_cg(ip)=( xgl(i)+1 )*dx/2 + x0;
      intma_cg(i,e)=ip;
   end %i
end %e

%Generate Periodicity pointer for CG
for i=1:npoin_cg
   periodicity_cg(i)=i;
end
periodicity_cg(npoin_cg)=periodicity_cg(1);

%Generate COORD and INTMA for DG
ip=0;
for e=1:nelem
    for i=1:ngl
        ip=ip+1;
        intma_dg(i,e)=ip;
    end
end

for e=1:nelem
   for i=1:ngl
      ip_cg=intma_cg(i,e);
      ip_dg=intma_dg(i,e);
      coord_dg(ip_dg)=coord_cg(ip_cg);
   end 
end

for i=1:npoin_dg
    periodicity_dg(i)=i;
end      