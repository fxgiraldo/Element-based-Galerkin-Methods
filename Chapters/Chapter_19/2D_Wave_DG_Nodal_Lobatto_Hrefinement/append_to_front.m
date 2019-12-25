%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [list,nlist] = append_to_front(lin1,nlin1,lin2,nlin2)

    list = [1:nlin1+nlin2];
    for i=1:nlin1
       list(i) = lin1(i);
    end
    
    for i=1:nlin2
       list(nlin1+i) = lin2(i); 
    end
    nlist = nlin1+nlin2;

end
