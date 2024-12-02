%---------------------------------------------------------------------%
%This function uses the Krividonova (Moment) Limiter
%Written by F.X. Giraldo on 1/2011
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [limit_element,qp,qmodal] = limiter_nodal_v2(qp,vdm,vdm_inv,intma,nelem,ngl)

%Initialize
limit_element=zeros(nelem,1);
qmodal=zeros(ngl,nelem);
qtemp=zeros(ngl,nelem);
umodal=zeros(ngl,nelem);
utemp=zeros(ngl,nelem);
    
for e=1:nelem
    for i=1:ngl
        ip=intma(i,e);
        qtemp(i,e)=qp(ip,1);
        utemp(i,e)=qp(ip,2);
    end
end

for e=1:nelem
   qmodal(:,e)=vdm_inv(:,:)*qtemp(:,e);
   umodal(:,e)=vdm_inv(:,:)*utemp(:,e);
end
[limit_element,qmodal] = limiter_modal(qmodal,nelem,ngl);
[limit_element,umodal] = limiter_modal(umodal,nelem,ngl);

for e=1:nelem
   qtemp(:,e)=vdm(:,:)*qmodal(:,e);
   utemp(:,e)=vdm(:,:)*umodal(:,e);
   for i=1:ngl
       ip=intma(i,e);
       qp(ip,1)=qtemp(i,e);
       qp(ip,2)=utemp(i,e);
   end
end
