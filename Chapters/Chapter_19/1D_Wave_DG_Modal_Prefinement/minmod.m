%---------------------------------------------------------------------%
%This is the MINMOD Functon
%Written by F.X. Giraldo on 1/2011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qlimit,limit] = minmod(q1,q2,q3)

tol=1e-16;

sign1=sign(q1);
sign2=sign(q2);
sign3=sign(q3);
q=zeros(3,1);
q(1)=abs(q1);
q(2)=abs(q2);
q(3)=abs(q3);

if (sign1 == sign2 && sign2 == sign3)
    qlimit=sign1*min(q);
else
    qlimit=0;
end

limit=1;
if ( abs(qlimit-q1) <= tol) 
    limit=0;
end
