%---------------------------------------------------------------------%
%This function computes the Mass Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = apply_dss(rhs,intma,Mmatrix_inv,npoin,nelem,ngl,iperiodic)

%Initialize
rhs_continuous=zeros(npoin,2);

%Gather Solution to make it Global
for e=1:nelem
    for i=1:ngl
        ip=iperiodic(intma(i,e));
        for k=1:2
            rhs_continuous(ip,k)=rhs_continuous(ip,k) + rhs(k,i,e);
        end
    end %i
end %e
   
%Scatter Solution to make it Local
rhs=0;
rhs_continuous(:,1)=Mmatrix_inv*rhs_continuous(:,1);
rhs_continuous(:,2)=Mmatrix_inv*rhs_continuous(:,2); 
for e=1:nelem
    for i=1:ngl
        ip=iperiodic(intma(i,e));
        for k=1:2
            rhs(k,i,e)=rhs_continuous(ip,k);
        end
    end %i
end %e   
      
