%---------------------------------------------------------------------%
%This function applies the Legendre Filter.
%Written by F.X. Giraldo on 7/2007
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = apply_filter2D(qp,f,intma,iperiodic,jac,wnq,Mmatrix,npoin,nelem,ngl)

%Initialize
rhs=zeros(npoin,1);
q_e=zeros(ngl,ngl);

%Integrate Divergence of Flux
for e=1:nelem
   
   %Store Values in Tensor-Product Form
   for j=1:ngl
       for i=1:ngl
           I=iperiodic(intma(e,i,j));
           q_e(i,j)=qp(I);
       end
   end
   q_f=f*q_e*f';
   
   %Store Values in Long-Vector Form
   for j=1:ngl
       for i=1:ngl
           wq=wnq(i)*wnq(j)*jac(e,i,j);
           I=iperiodic(intma(e,i,j));
           rhs(I)=rhs(I) + wq*q_f(i,j);
       end
   end

end %e

%Multiply by Inverse Mass
rhs=rhs./Mmatrix;

%Periodicity
for i=1:npoin
    rhs(i)=rhs(iperiodic(i));
end
