%---------------------------------------------------------------------%
%this function renumbers points on the edge between two points 
% e1, e2 are ngl x ngl intma entries
% f gives the local face numbers
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [e1] = renumber_face_points(e1,e2,f2)

ngl = size(e1,2);

switch f2
    case 1
        temp = e2(1,1:ngl,1);
        f1 = 2;
    case 2
        temp = e2(1,1:ngl,ngl);
        f1 = 1;
    case 3
        temp = e2(1,1,1:ngl);
        f1 = 4;
    case 4
        temp = e2(1,ngl,1:ngl);
        f1 = 3;
end

switch f1
    case 1
        e1(1,1:ngl,1) = temp;
    case 2
        e1(1,1:ngl,ngl) = temp;
    case 3
        e1(1,1,1:ngl) = temp;
    case 4
        e1(1,ngl,1:ngl) = temp;
end
