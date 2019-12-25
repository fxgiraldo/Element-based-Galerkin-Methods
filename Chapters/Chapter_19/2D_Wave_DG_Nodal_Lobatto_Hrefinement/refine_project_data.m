%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [q] = refine_project_data(q,iref,nref,P1s2d,P2s2d,P3s2d,P4s2d,ngl,tc)

for i=1:nref
    
    qa = zeros(1,ngl*ngl);   
    qb = zeros(1,ngl*ngl);


     for j=1:ngl
        for k=1:ngl
            n=(j-1)*ngl+k;
            qa(n) = q(iref(i),k,j);
        end
    end
    
    qb = qa*P1s2d;
    for j=1:ngl
        for k=1:ngl
            n=(j-1)*ngl+k;
            q(tc(1,iref(i)),k,j) = qb(n);
        end
    end
    
    qb = qa*P2s2d;
    for j=1:ngl
        for k=1:ngl
            n=(j-1)*ngl+k;
            q(tc(2,iref(i)),k,j) = qb(n);
        end
    end
    
    qb = qa*P3s2d;
    for j=1:ngl
        for k=1:ngl
            n=(j-1)*ngl+k;
            q(tc(3,iref(i)),k,j) = qb(n);
        end
    end
    
    qb = qa*P4s2d;
    for j=1:ngl
        for k=1:ngl
            n=(j-1)*ngl+k;
            q(tc(4,iref(i)),k,j) = qb(n);
        end
    end
    
end
