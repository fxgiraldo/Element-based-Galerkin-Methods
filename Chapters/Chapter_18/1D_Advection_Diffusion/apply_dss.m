%---------------------------------------------------------------------%
%This function computes the Mass Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = apply_dss(rhs,intma,Mmatrix_inv,npoin,nelem,ngl,iperiodic)

%Initialize
rhs_continuous=zeros(npoin,1);

%Gather Solution to make it Global
for e=1:nelem
    for i=1:ngl
        ip=iperiodic(intma(i,e));
        rhs_continuous(ip)=rhs_continuous(ip) + rhs(i,e);
    end %i
end %e
   
%Scatter Solution to make it Local
%rhs=0;
rhs_continuous=Mmatrix_inv*rhs_continuous;
for e=1:nelem
    for i=1:ngl
        ip=iperiodic(intma(i,e));
        rhs(i,e)=rhs_continuous(ip);
    end %i
end %e   
      
