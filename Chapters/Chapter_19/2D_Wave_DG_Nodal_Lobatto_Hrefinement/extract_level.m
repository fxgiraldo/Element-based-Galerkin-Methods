%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [list,nlist] = extract_level(lin,nlin,tl,lev)

list = 0;
k=0;
for i=1:nlin

   if tl(lin(i))==lev
        k=k+1;
        list(k) = lin(i);
   end
end
nlist = k;

end
