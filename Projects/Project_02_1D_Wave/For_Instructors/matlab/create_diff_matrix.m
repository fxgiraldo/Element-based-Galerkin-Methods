%---------------------------------------------------------------------%
%This function computes the Differentiation matrix.
%Written by F.X. Giraldo on 5/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function element_diff = create_diff_matrix(ngl,nq,wnq,psi,dpsi)

%Initialize
element_diff=zeros(ngl,ngl);
    
%LGL Integration
for k=1:nq
  wk=wnq(k);
  for i=1:ngl
      dhdx_i=dpsi(i,k);
      for j=1:ngl
          h_j=psi(j,k);
          element_diff(i,j)=element_diff(i,j) + wk*dhdx_i*h_j;
      end %j
  end %i
end %k      
