%---------------------------------------------------------------------%
%This function applies the DSS operation (Gather then a Scatter)
%Written by F.X. Giraldo on October 26, 2011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [rhs] = apply_dss(rhs,intma,iperiodic,ngl,npoin,nelem)

%Initialize array
rhs_continuous=zeros(npoin,1);

%Gather Solution (Q^T operator) to make it Global
for e=1:nelem
 for j=1:ngl
    for i=1:ngl
       ip=iperiodic( intma(e,i,j) );
       rhs_continuous(ip)=rhs_continuous(ip) + rhs(e,i,j);
    end %i
 end %k
end %e

%Scatter Solution (Q operator) to make it Local
rhs=0;
for e=1:nelem
 for j=1:ngl
    for i=1:ngl
       ip=iperiodic( intma(e,i,j) );
       rhs(e,i,j)=rhs_continuous(ip);
    end %i
 end %k
end %e
 


      
