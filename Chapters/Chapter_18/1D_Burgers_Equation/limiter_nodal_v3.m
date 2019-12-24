%---------------------------------------------------------------------%
%This function uses the Krividonova (Moment) Limiter
%Written by F.X. Giraldo on 5/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [limit_element,qp,qmodal] = limiter_nodal_v3(qp,vdm,vdm_inv,intma,nelem,ngl,limit)

%Initialize
limit_element=zeros(nelem,1);
qmodal=zeros(ngl,nelem);
qtemp=zeros(ngl,nelem);
qnodal=zeros(ngl,1);

for k=2:-1:1
    
    %Store Element-wise Solution
    for e=1:nelem
        for i=1:ngl
            ip=intma(i,e);
            qtemp(i,e)=qp(ip,k);
        end
    end
    
    %Compute Modal Solution
    for e=1:nelem
       qmodal(:,e)=vdm_inv(:,:)*qtemp(:,e);
    end
    
    %Limit Solution
    if (limit > 0)
        [limit_element,qmodal] = limiter_modal(qmodal,nelem,ngl);
    end
    
    %Compute Nodal Solution
    for e=1:nelem
        qnodal(:)=vdm(:,:)*qmodal(:,e);
        for i=1:ngl
            ip=intma(i,e);
            qp(ip,k)=qnodal(i);
        end
    end
end %k

%Get Modal Solution for Plotting
for e=1:nelem
    for i=1:ngl
        ip=intma(i,e);
        qtemp(i,e)=qp(ip,1);
    end
end
for e=1:nelem
   qmodal(:,e)=vdm_inv(:,:)*qtemp(:,e);
end
