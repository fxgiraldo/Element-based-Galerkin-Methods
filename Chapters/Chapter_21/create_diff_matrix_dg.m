%---------------------------------------------------------------------%
%This function computes the Differentiation matrix.
%Written by F.X. Giraldo on 5/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function Dmatrix = create_diff_matrix_dg(ngl,nq,wnq,psi,dpsi)

%Initialize
Dmatrix=zeros(ngl,ngl);

%Integrate Divergence of Flux
%for ie=1:nelem
    
   %LGL Integration
   for i=1:ngl
       for j=1:ngl
           for k=1:nq
               wk=wnq(k);
               dhdx_ik=dpsi(i,k);
               h_jk=psi(j,k);
               Dmatrix(i,j)=Dmatrix(i,j) + wk*dhdx_ik*h_jk;
           end %k
       end %j
   end %i
%end %ie

      
