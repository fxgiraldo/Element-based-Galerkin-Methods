%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [list] = level_sort(lin,nlin,tl,maxlev,type)
k=0;
if type==1
    for l=1:maxlev

        for i=1:nlin

           if tl(lin(i))==l
                k=k+1;
                list(k) = lin(i);
           end
        end
    end
elseif type==2
    for l=maxlev:-1:0

        for i=1:nlin

           if tl(lin(i))==l
                k=k+1;
                list(k) = lin(i);
           end
        end
    end
end
end
