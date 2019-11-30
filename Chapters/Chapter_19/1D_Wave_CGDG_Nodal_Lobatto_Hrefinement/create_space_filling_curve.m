%---------------------------------------------------------------------%
%This function creates the space-filling curve
%Written by F.X. Giraldo on 3/13/2016 using MA Kopera's matlab algorithm
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [sfc,nsfc] = create_space_filling_curve(nelem,nelem0,children,active)

%initialize
sfc=zeros(nelem,1);
m=0;

for e=1:nelem0
    [ptr,sfc,m] = is_in_tree(sfc,children,active,e,m);
end
nsfc=m;