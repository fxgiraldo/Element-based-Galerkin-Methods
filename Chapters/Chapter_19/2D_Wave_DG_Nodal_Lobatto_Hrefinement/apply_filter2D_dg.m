%---------------------------------------------------------------------%
%This function applies the Legendre Filter.
%Written by F.X. Giraldo on 7/2007
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = apply_filter2D_dg(qp,f,nelemT,sfc,nsfc,ngl)

%Initialize
rhs=zeros(nelemT,ngl,ngl);
q_e=zeros(ngl,ngl);

%Integrate Divergence of Flux
for is=1:nsfc
    e=sfc(is);
   
   %Store Coordinates
    q_e(:,:)=qp(e,:,:);
    rhs(e,:,:)=f*q_e*f';

end %e

