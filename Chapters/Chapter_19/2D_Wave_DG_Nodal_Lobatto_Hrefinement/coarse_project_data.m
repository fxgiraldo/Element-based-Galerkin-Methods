%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [q] = coarse_project_data(q,icor,ncor,P1g2d,P2g2d,P3g2d,P4g2d,ngl,tc,tp)

for i=1:ncor
    
    qa = zeros(1,ngl*ngl);
    qb = zeros(1,ngl*ngl);
    qc = zeros(1,ngl*ngl);
    qd = zeros(1,ngl*ngl);

    ch1 = tc(1,tp(icor(i)));
    ch2 = tc(2,tp(icor(i)));
    ch3 = tc(3,tp(icor(i)));
    ch4 = tc(4,tp(icor(i)));
    
    for j=1:ngl
        for k=1:ngl
            n=(j-1)*ngl+k;
            
                qa(n) = q(ch1,k,j);   
                qb(n) = q(ch2,k,j);    
                qc(n) = q(ch3,k,j); 
                qd(n) = q(ch4,k,j);
        end
    end


    qa = qa*P1g2d + qb*P2g2d + qc*P3g2d + qd*P4g2d;
    
    for j=1:ngl
        for k=1:ngl
            n=(j-1)*ngl+k;
            q(tp(icor(i)),k,j) = qa(n);
        end


    end
end
end
