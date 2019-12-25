%---------------------------------------------------------------------%
%This function applies the DSS operation (Gather then a Scatter)
%Written by F.X. Giraldo on October 26, 2011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [rhs] = apply_dss(rhs,intma,iperiodic,mass_global,ngl,npoinT,sfc,nsfc,face,ffc,nffc)

%Initialize array
rhs_continuous=zeros(npoinT,1);

%gather from hanging nodes


%Gather Solution (Q^T operator) to make it Global
for is=1:nsfc
    e=sfc(is);
 for j=1:ngl
    for i=1:ngl
       ip=iperiodic( intma(e,i,j) );
       rhs_continuous(ip)=rhs_continuous(ip) + rhs(e,i,j);
    end %i
 end %k
end %e

rhs_continuous = gather_from_hanging_nodes_CG_DG(rhs,rhs_continuous,intma,face,ffc,nffc,ngl);



%Solve globally and scatter Solution (Q operator) to make it Local
rhs=0;
for is=1:nsfc
    e=sfc(is);
 for j=1:ngl
    for i=1:ngl
       ip=iperiodic( intma(e,i,j) );
       rhs(e,i,j)=rhs_continuous(ip)/mass_global(ip);
    end %i
 end %k
end %e

rhs = scatter_to_hanging_nodes_CG_DG(rhs,face,ffc,nffc,ngl);

end
 


      
