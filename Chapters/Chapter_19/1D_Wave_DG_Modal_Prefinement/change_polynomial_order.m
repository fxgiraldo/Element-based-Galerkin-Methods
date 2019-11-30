%---------------------------------------------------------------------%
%This function uses the Krividonova (Moment) Limiter
%Written by F.X. Giraldo on 1/2011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qp,element_order] = change_polynomial_order(qp,nelem,element_order,p_tol_c,p_tol_r,p_tol2,nop_min,nop_max)

%P-Derefinement
for e=1:nelem
   N=element_order(e)+1;
   i=N;
   denominator=norm(qp(1:N,e));
   limit=1;
   
   while limit == 1

       %Get Coefficients
       numerator=norm(qp(i:N,e));
             
       %Check if current mode is significant
       check=numerator/denominator;
       if ( check <= p_tol_c || denominator <= p_tol2)
            if i == nop_min+1
                limit = 0; %stop because reached limit
            else
                i=i-1; %go to the next lower mode
                if (i < nop_min+1)
                    i=nop_min+1;
                    limit = 0;
                end
            end
            disp(['(1) Derefine: e= ',num2str(e),' p= ',num2str(i),' num = ',num2str(numerator),' den= ',num2str(denominator),' check= ',num2str(check)]);
       else
           limit = 0; %stop because mode is significant
       end %if
       
   end %while
   element_order(e)=i-1;
end %e

%P-Enrichment
for e=1:nelem
   N=element_order(e)+1;
   i=N;
   denominator=norm(qp(1:N,e));
   
   %Get Coefficients
   numerator=norm(qp(i:N,e));

   %Check if current mode is significant
   check=numerator/denominator;
   if (check > p_tol_r && denominator > p_tol2)
       if (i<nop_max+1)
           i=i+1; %need to enrich the space
           disp(['(2) Enrich: e= ',num2str(e),' p= ',num2str(i),' num = ',num2str(numerator),' den= ',num2str(denominator),' check= ',num2str(check)]);
       end %if
   end %if   
   element_order(e)=i-1;
end %e
