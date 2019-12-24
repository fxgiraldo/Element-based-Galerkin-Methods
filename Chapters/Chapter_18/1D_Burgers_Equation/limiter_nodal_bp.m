%---------------------------------------------------------------------%
%This function uses the Krividonova (Moment) Limiter
%Written by F.X. Giraldo on 5/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [limit_element,qp,qmodal] = limiter_nodal_bp(qp,q0,vdm,vdm_inv,intma,nelem,ngl)

%Initialize
limit_element=zeros(nelem,1);
qtemp=zeros(ngl,nelem);
qnodal=zeros(ngl,nelem);
qmodal=zeros(ngl,nelem);
qmodal_loc=zeros(ngl,1);
qlocal=zeros(ngl,1);
eps=1e-6;
tvector=zeros(3,1);

%Loop through variables
for k=1:2
    qtemp=zeros(ngl,nelem);
    for e=1:nelem
        for i=1:ngl
            ip=intma(i,e);
            qtemp(i,e)=q0(ip,k);
            qnodal(i,e)=qp(ip,k);
        end
    end

    %Global Extrema
    global_max=max(max(qtemp));
    global_min=min(min(qtemp));

    %Loop through element
    for e=1:nelem

        %Store Local Space
        qlocal(:)=qnodal(:,e);

        %Map to Modal Space
        qmodal_loc(:)=vdm_inv(:,:)*qlocal(:);

        %Local Extrema
        local_max=max(qlocal);
        local_min=min(qlocal);

        %Computer Limiter Weight
        qmean=qmodal_loc(1);
        t1_den=local_max - qmean;
        t2_den=local_min - qmean;
        tvector(1)=abs((global_max - qmean)/( (local_max - qmean) + eps ));
        tvector(2)=abs((global_min - qmean)/( (local_min - qmean) + eps ));
        tvector(3)=1;
        theta=min(tvector);
        if (theta < 1)
            limit_element(e)=1;
        end

        for i=1:ngl
            ip=intma(i,e);
            qp(ip,k)=theta*qlocal(i) + (1-theta)*qmean;
        end
    end %e
end %k

%Map to Modal Space
qtemp=zeros(ngl,nelem);
for e=1:nelem
    for i=1:ngl
        ip=intma(i,e);
        qtemp(i,e)=qp(ip,1);
    end
end

for e=1:nelem
   qmodal(:,e)=vdm_inv(:,:)*qtemp(:,e);
end



