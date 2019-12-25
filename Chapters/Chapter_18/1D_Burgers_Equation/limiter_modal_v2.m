%---------------------------------------------------------------------%
%This function uses the Krividonova (Moment) Limiter
%Written by F.X. Giraldo on 5/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [limit_element,qp] = limiter_modal(qp,nelem,ngl)

%Initialize
limit_element=zeros(nelem,1);

%Integrate Flux Terms
for ie=1:nelem
   iel=ie-1;
   ier=ie+1;
   if (ie == 1) 
      iel=nelem;
   end 
   if (ie == nelem) 
      ier=1;
   end 
   iend=ngl;
   iterations=0;
   limit=1;
   
   while limit == 1
       
       iterations=iterations + 1;

       %Get Coefficients
       q_e=qp(iend,ie);
       dq_p=qp(iend-1,ier) - qp(iend-1,ie);
       dq_m=qp(iend-1,ie)  - qp(iend-1,iel);
       
       %Limit Solution
       [qlimit,limit]=minmod(q_e,dq_p,dq_m);
       
       if iterations == 1
           limit_element(ie)=limit;
       end
       
       %Limit Solution
       if limit == 1 
            qp(iend,ie)= qlimit;
            iend=iend - 1;
            %iend_min=round(ngl/2);
            iend_min=1;
            if iend == iend_min
                limit = 0;
            end
       end %if
       
   end %while
   
end %ie
