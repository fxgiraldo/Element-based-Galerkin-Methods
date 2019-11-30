%---------------------------------------------------------------------%
%This function uses the Krividonova (Moment) Limiter
%Written by F.X. Giraldo on 1/2011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [limit_element,qp] = limiter_modal_krividonova(qp,nelem,nq,element_order)

%Initialize
limit_element=zeros(nelem,1);

%Integrate Flux Terms
for e=1:nelem
   el=e-1;
   er=e+1;
   if (e == 1) 
      el=nelem;
   end 
   if (e == nelem) 
      er=1;
   end 
   Ne=element_order(e)+1;
   Nl=element_order(el)+1;
   Nr=element_order(er)+1;
   iend_e=Ne;
   iend_l=Nl;
   iend_r=Nr;
   iterations=0;
   limit=1;
   
   while limit == 1
       
       iterations=iterations + 1;

       %Get Coefficients
       q_e=qp(iend_e,e);
       dq_p=qp(iend_r-1,er) - qp(iend_e-1,e);
       dq_m=qp(iend_e-1,e)  - qp(iend_l-1,el);
       
       %Limit Solution
       [qlimit,limit]=minmod(q_e,dq_p,dq_m);
       
       if iterations == 1
           limit_element(e)=limit;
       end
       
       %Limit IEND Mode of the Solution
       if limit == 1 
            qp(iend_e,e)= qlimit;
            iend_e=iend_e - 1;
            %iend_min=round(ngl/2);
            iend_min=1;
            if iend == iend_min
                limit = 0;
            end
       end %if
       
   end %while
   
end %e
