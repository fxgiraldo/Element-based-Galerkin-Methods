%---------------------------------------------------------------------%
%This function applies the DSS operation (Gather then a Scatter)
%Written by F.X. Giraldo on October 26, 2011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [q] = apply_dss_vector(q,DG_to_CG,npoin,npoin_CG)

%Initialize array
q_continuous=zeros(npoin_CG,1);

%Gather Solution (Q^T operator) to make it Global
for i=1:npoin
    ip_CG=DG_to_CG(i);
    q_continuous(ip_CG)=q_continuous(ip_CG) + q(i);      
end %i

%Scatter Solution (Q operator) to make it Local
q=zeros(npoin,1);
for i=1:npoin
    ip_CG=DG_to_CG(i); 
    q(i)=q_continuous(ip_CG);
end %i
% q_continuous