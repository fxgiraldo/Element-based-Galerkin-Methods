%---------------------------------------------------------------------%
%This code computes the minmod funtion
%Written by F.X. Giraldo on 4/2000
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function qp_mod = minmod(qp,delta_plus,delta_minus,dx,m)

%Constants
%m=0.35;
c=m*dx*dx;
if delta_plus == 0
    delta_plus=-1e-10;
end
if delta_minus == 0
    delta_minus=-1e-10;
end
if ( abs(qp) <= abs(c))
    
   qp_mod=qp;
   
   elseif ( (sign(qp) == sign(delta_plus) ))%== sign(delta_minus)))% && (sign(qp) == sign(delta_minus)) )
    if(sign(qp) == sign(delta_minus))
    t=zeros(1,3);
    t(1)=abs(qp);
    t(2)=abs(delta_plus);
    t(3)=abs(delta_minus);
    qp_mod=sign(qp)*min(t)*0.5; %works only for linear elements now
    else
        qp_mod=0;
    end
   
else
    qp_mod=0;
    
end
      
