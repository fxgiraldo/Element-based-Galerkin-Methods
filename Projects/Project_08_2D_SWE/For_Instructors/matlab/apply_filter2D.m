%---------------------------------------------------------------------%
%This function applies the Legendre Filter.
%Written by F.X. Giraldo on 7/2007
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = apply_filter2D(qp,f,intma,iperiodic,jac,wnq,Mmatrix,npoin,nelem,ngl)

%Initialize
rhs=zeros(3,npoin);
q_e=zeros(3,ngl,ngl);
q_f=zeros(3,ngl,ngl);
qq=zeros(ngl,ngl);

%Integrate Divergence of Flux
for e=1:nelem
   
   %Store Values in Tensor-Product Form
   for j=1:ngl
       for i=1:ngl
           I=iperiodic(intma(i,j,e));
           for m=1:3
            q_e(m,i,j)=qp(m,I);
           end
       end
   end
   for m=1:3
       qq(:,:)=q_e(m,:,:);
       q_prod=f*qq*f';
       q_f(m,:,:)=q_prod(:,:);
   end

   %Store Values in Long-Vector Form
   for j=1:ngl
       for i=1:ngl
           wq=wnq(i)*wnq(j)*jac(i,j,e);
           I=iperiodic(intma(i,j,e));
           for m=1:3
               rhs(m,I)=rhs(m,I) + wq*q_f(m,i,j);
           end
       end
   end

end %e

%Multiply by Inverse Mass
for i=1:npoin
    for m=1:3
        rhs(m,i)=rhs(m,i)/Mmatrix(i);
    end
end

%Periodicity
for i=1:npoin
    for m=1:3
        rhs(m,i)=rhs(m,iperiodic(i));
    end
end
