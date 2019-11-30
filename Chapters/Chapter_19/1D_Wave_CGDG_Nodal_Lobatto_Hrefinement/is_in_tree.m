%---------------------------------------------------------------------%
%This function finds, recursively, the child of a root element that is
%active.
%Modified by F.X. Giraldo on 3/13/2016 to use Binary-trees from Michal Kopera's Quad-tree routine
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [ptr,sfc,m] = is_in_tree(sfc,children,active,ptr,m)

if (active(ptr)==1)
    m=m+1;
    sfc(m) = ptr;
else
    if(children(1,ptr)~=0)
       for i=1:2
          [ptr1,sfc,m] = is_in_tree(sfc,children,active,children(i,ptr),m);
       end
    end
end