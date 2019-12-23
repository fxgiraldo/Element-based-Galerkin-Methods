%---------------------------------------------------------------------%
%This function applies the DSS operation (Gather then a Scatter)
%Written by F.X. Giraldo on October 26, 2011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = apply_dss_rhs(rhs,DG_to_CG,npoin,npoin_CG)

%Initialize array
rhs_continuous=zeros(npoin_CG,1);

%Gather Solution (Q^T operator) to make it Global
for i=1:npoin
    ip_CG=DG_to_CG(i);
    rhs_continuous(ip_CG)=rhs_continuous(ip_CG) + rhs(i);
end %i

%Scatter Solution (Q operator) to make it Local
rhs=zeros(npoin,1);
icounter_CG=zeros(npoin_CG,1);
for i=1:npoin
    ip_CG=DG_to_CG(i); 
    icounter_CG(ip_CG)=icounter_CG(ip_CG) + 1;
    if icounter_CG(ip_CG) == 1
        rhs(i)=rhs_continuous(ip_CG);
    elseif icounter_CG(ip_CG) > 1
        rhs(i)=0;
    end
end %i