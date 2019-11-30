%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [Me, Dwe, Fe] = element_matrices(psi,dpsi,ngl,nq,wnq)

%Initialize
Me=zeros(ngl,ngl);
Dwe=zeros(ngl,ngl);
Fe=zeros(ngl,ngl);
 
%Do LGL Integration
for l=1:nq
  wq=wnq(l);
  for j=1:ngl
     for i=1:ngl
        Me(i,j)=Me(i,j) + wq*psi(i,l)*psi(j,l);
        Dwe(i,j)=Dwe(i,j) + wq*dpsi(i,l)*psi(j,l);
     end %j
  end %i
end %l
Fe(1,1)=-1;
Fe(ngl,ngl)=+1;



      
