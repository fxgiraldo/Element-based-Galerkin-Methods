%---------------------------------------------------------------------%
%This function uses the Krividonova (Moment) Limiter
%Written by F.X. Giraldo on 1/2011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qp,element_order] = change_polynomial_order(qp,VDM_Matrix,VDM_Inv_Matrix,nelem,element_order,p_tol_c,p_tol_r,p_tol2,nop_min,nop_max)
eps=1e-15;
qm=zeros(nop_max+1,1);

%P-Adapt
for e=1:nelem
   N=element_order(e)+1;
   i=N;   
   iadapt=0;
   
   %Get Modal Coefficients
   qm(1:nop_max+1)=0;
   qm(1:i)=VDM_Inv_Matrix(1:i,1:i,i-1)*qp(1:i,e);

   %Check Criterion
   numerator=norm(qm(i:N));
   denominator=norm(qm(1:N))+eps;
   check=numerator/denominator;
   
   %Check if current mode is significant
   if (check <= p_tol_c || denominator <= p_tol2) %Coarsen
       if (i>=nop_min+2)
           i=i-1; %need to coarsen the space
           iadapt=1;
%            disp(['(1) Coarsen: e= ',num2str(e),' p= ',num2str(i),' num = ',num2str(numerator),' den= ',num2str(denominator),' check= ',num2str(check)]);
       end %if
       element_order(e)=i-1;
       %Get Nodal Coefficients
       qp(1:i,e)=VDM_Matrix(1:i,1:i,i-1)*qm(1:i);
   elseif (check > p_tol_r && denominator > p_tol2) %Refine
       if (i<nop_max+1)
           i=i+1; %need to enrich the space
           iadapt=1;
%            disp(['(2) Enrich: e= ',num2str(e),' p= ',num2str(i),' num = ',num2str(numerator),' den= ',num2str(denominator),' check= ',num2str(check)]);
       end %if
       element_order(e)=i-1;
       %Get Nodal Coefficients
       qp(1:i,e)=VDM_Matrix(1:i,1:i,i-1)*qm(1:i);
   end %if      

   %Special condition for Linear
   if (iadapt == 0 && i==2 && denominator > p_tol2) %Refine
       if (i<nop_max+1)
           i=i+1; %need to enrich the space
%            disp(['(2) Enrich: e= ',num2str(e),' p= ',num2str(i),' num = ',num2str(numerator),' den= ',num2str(denominator),' check= ',num2str(check)]);
       end %if
       element_order(e)=i-1;
       %Get Nodal Coefficients
       qp(1:i,e)=VDM_Matrix(1:i,1:i,i-1)*qm(1:i);
   end %if      
end %e
