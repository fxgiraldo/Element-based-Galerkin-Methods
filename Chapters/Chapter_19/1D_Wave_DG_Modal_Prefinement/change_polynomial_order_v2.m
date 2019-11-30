%---------------------------------------------------------------------%
%This function uses the Krividonova (Moment) Limiter
%Written by F.X. Giraldo on 1/2011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qp,element_order] = change_polynomial_order_v2(qp,nelem,element_order,p_tol_c,p_tol_r,p_tol2,nop_min,nop_max)
eps=1e-15;

%P-Adapt
for e=1:nelem
   N=element_order(e)+1;
   i=N;
   denominator=norm(qp(1:N,e))+eps; %Norm of All modes
   
   %Get Coefficients
   numerator=norm(qp(i:N,e)); %Norm of Nth mode

   %Check if current mode is significant
   check=numerator/denominator;
   if (check <= p_tol_c || denominator <= p_tol2)
       if (i>=nop_min+2)
           i=i-1; %need to coarsen the space
           disp(['(1) Coarsen: e= ',num2str(e),' p= ',num2str(i),' num = ',num2str(numerator),' den= ',num2str(denominator),' check= ',num2str(check)]);
       end %if
   elseif (check > p_tol_r && denominator > p_tol2)
       if (i<nop_max+1)
          i=i+1; %need to enrich the space
          disp(['(2) Enrich: e= ',num2str(e),' p= ',num2str(i),' num = ',num2str(numerator),' den= ',num2str(denominator),' check= ',num2str(check)]);
       end %if
   end %if   
   element_order(e)=i-1;
end %e
