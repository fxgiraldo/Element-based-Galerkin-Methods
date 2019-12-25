%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [icorc,ncorc] = correct_coarse_list_new(icor,ncor,tp,tc,tm)

%sort
icor = sort(icor);
icorc = 0;
ncorc = 0;
i=1;
m=1;
k=0;
parent_lock = 0;
if ncor >=4
while i<=ncor-3

    if icor(i)==tc(1,tp(icor(i)))
        parent_lock = tp(icor(i));
        m=m+1;
        for j=1:3
            if icor(i+j)==tc(m,parent_lock) 
                m=m+1;
            else
                m=1;
                i=i+j;
                break
            end
        end
        
        if m==5
%             if tm(icor(i))==1
                k=k+1;
                icorc(k) = icor(i);
%             end
            m=1;
            i=i+4;
        end
    else
        i=i+1;
        
    end
    
%     if icor(i)==tc(m,tp(icor(i))) 
%         i=i+1;
%         m=m+1;
%     else
% 
%         m=1;
%         i=i+1;
%         continue
%     end
%     
%     if m==5
%         k=k+1;
%         icorc(k) = icor(i-1);
%         m=1;
%     end
end


ncorc = k;
end
end
