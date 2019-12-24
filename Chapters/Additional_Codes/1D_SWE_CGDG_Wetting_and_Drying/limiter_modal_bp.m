%---------------------------------------------------------------------%
%This function uses the Bound-Preserving Limiter
%Written by F.X. Giraldo on 3/2013
%           Department of Applied Maths
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function qp = limiter_modal_bp(qp,qb,vdm,vdm_inv,nelem,ngl)

%Initialize
limit_element=zeros(nelem,1);
qe_ave=zeros(nelem,1);
qe_modal=zeros(ngl,1);
weight=zeros(ngl,1);
qtemp=zeros(ngl,1);

%Compute Average/Mean Values of all elements
for e=1:nelem
   temp=0;
   for i=1:ngl
    qtemp(i)=qp(1,i,e) + qb(i,e);
    temp=temp + qtemp(i);
   end
   %qe_ave(e)=vdm_inv(1,:)*qtemp;
   qe_ave(e)=temp/ngl;
end %e

%Compute Global Extrema
% global_min=min(qp(:));
% global_max=max(qp(:));
global_min=0;
global_max=0.6;

%Loop through All Elements 
for e=1:nelem

    for i=1:ngl
     qtemp(i)=qp(1,i,e) + qb(i,e);
    end
   
    %Compute local extrema
   local_min=min(qtemp(:));
   local_max=max(qtemp(:));
   
   %compute linear coefficient/weight
    local_mean=qe_ave(e);
    t1=(global_max - local_mean)/(local_max - local_mean);
    t2=(global_min - local_mean)/(local_min - local_mean);
   alpha=1;
   alpha=min(alpha,t1);
   alpha=min(alpha,t2);
   
   for i=1:ngl
    qp(1,i,e)=alpha*( qtemp(i) - local_mean ) + local_mean - qb(i,e);
   end
   
end %e
