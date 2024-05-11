%---------------------------------------------------------------------%
%This function computes the Mass Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = apply_dss(rhs,intma,intma_cg,coord,Mmatrix,npoin,npoin_cg,nelem,ngl,nq,wnq,psi,iperiodic,ivolume_integrate)

%Initialize
rhs_continuous=zeros(npoin_cg,2);

%Volume Integrate
if (ivolume_integrate == 1)
    rhs_volume = compute_volume_integral(rhs,intma,coord,npoin,nelem,ngl,nq,wnq,psi);
else
    rhs_volume=rhs;
end

%Gather Solution to make it Global
for e=1:nelem
    for i=1:ngl
        ip_cg=intma_cg(i,e);
        ip=iperiodic(intma(i,e));
        for k=1:2
            rhs_continuous(ip_cg,k)=rhs_continuous(ip_cg,k) + rhs_volume(ip,k);
        end
    end %i
end %e
   
%Scatter Solution to make it Local
rhs=zeros(npoin,2);
rhs_continuous(:,1)=Mmatrix\rhs_continuous(:,1);
rhs_continuous(:,2)=Mmatrix\rhs_continuous(:,2); 
for e=1:nelem
    for i=1:ngl
        ip_cg=intma_cg(i,e);
        ip=iperiodic(intma(i,e));
        for k=1:2
            rhs(ip,k)=rhs_continuous(ip_cg,k);
        end
    end %i
end %e   
      
