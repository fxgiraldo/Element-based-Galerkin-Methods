%---------------------------------------------------------------------%
%This function uses a Zeroeth Order (mean value) FVM-type Limiter
%Written by F.X. Giraldo on 1/2011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [limit_element,qp] = limiter_modal_fvm(qp,nelem,ngl)

%Initialize
limit_element=zeros(nelem,1);
q_ave=zeros(nelem,1);
q_nodal=zeros(ngl,1);

%Store Mean Value
q_ave(:)=qp(1,:);

%Loop through left and right neighbors of Elements

for ie=1:nelem
   iel=ie-1;
   ier=ie+1;
   if (ie == 1) 
      iel=nelem;
   end 
   if (ie == nelem) 
      ier=1;
   end 
   
   plmin=q
   
end %ie
