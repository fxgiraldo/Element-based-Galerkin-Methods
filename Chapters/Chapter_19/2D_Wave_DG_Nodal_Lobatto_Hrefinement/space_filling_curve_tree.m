%---------------------------------------------------------------------%
%This subroutine creates the space filling curve out of the element tree
% note that nelem should be the original number of elements at the top
% level
% elmptr is the element pointer array which defines the element tree
% recursively
% nsfc is the length of sfc
%Written by M.A. Kopera on 10/2011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%

function [sfc,nsfc] = space_filling_curve_tree(nelem,tc,tm)

clear sfc;
sfc = 0;
m=0;
for e=1:nelem
[ptr,sfc,m] = is_in_tree(sfc,tm,tc,e,m);
end
nsfc=m;

end