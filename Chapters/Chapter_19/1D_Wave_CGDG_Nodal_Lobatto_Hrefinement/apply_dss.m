%---------------------------------------------------------------------%
%This function computes the Mass Matrix.
%Written by F.X. Giraldo on May 1, 2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = apply_dss(rhs,intma,npoin,nelem,element_order,active)

%Initialize
rhs_continuous=zeros(npoin,1);

%Gather Solution to make it Global
for e=1:nelem
    if (active(e) == 1)
        i=element_order(e);
        ii=i+1;
        for i=1:ii
            ip=intma(i,e);
            if (e == nelem && i == ii)
                ip = 1;
            end
            rhs_continuous(ip)=rhs_continuous(ip) + rhs(i,e);
        end %i
    end %active
end %e
   
%Scatter Solution to make it Local
for e=1:nelem
    if (active(e) == 1)
        i=element_order(e);
        ii=i+1;
        for i=1:ii
            ip=intma(i,e);
            if (e == nelem && i == ii)
                ip = 1;
            end
            rhs(i,e)=rhs_continuous(ip);
        end %i
    end %active
end %e   
      
