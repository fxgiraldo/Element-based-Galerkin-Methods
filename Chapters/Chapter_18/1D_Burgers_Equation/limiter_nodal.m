%---------------------------------------------------------------------%
%This function uses the Krividonova (Moment) Limiter
%Written by F.X. Giraldo on 1/2011
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [limit_element,qp] = limiter_nodal(qp,vdm,vdm_inv,nelem,ngl)

%Initialize
limit_element=zeros(nelem,1);
qmodal=zeros(ngl,nelem);
qtemp=zeros(ngl,nelem);

for k=1:2
    
    qtemp(:,:)=qp(k,:,:);

    for e=1:nelem
       qmodal(:,e)=vdm_inv(:,:)*qtemp(:,e);
    end
    [limit_element,qmodal] = limiter_modal(qmodal,nelem,ngl);
    for e=1:nelem
       qp(k,:,e)=vdm(:,:)*qmodal(:,e);
    end

end %k