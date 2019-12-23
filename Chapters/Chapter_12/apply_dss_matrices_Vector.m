%---------------------------------------------------------------------%
%This function applies the DSS operation (Gather then a Scatter)
%Written by F.X. Giraldo on October 26, 2011
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [Lmatrix,Rvector] = apply_dss_matrices_Vector(Lmatrix,Rvector,DG_to_CG,npoin,npoin_CG)

%Initialize array
Lmatrix_continuous=zeros(npoin_CG,npoin_CG);
Rvector_continuous=zeros(npoin_CG,1);

%Gather Solution (Q^T operator) to make it Global
for i=1:npoin
    ip_CG=DG_to_CG(i);
    Rvector_continuous(ip_CG)=Rvector_continuous(ip_CG) + Rvector(i);
    for j=1:npoin
        jp_CG=DG_to_CG(j);
        Lmatrix_continuous(ip_CG,jp_CG)=Lmatrix_continuous(ip_CG,jp_CG) + Lmatrix(i,j);
    end
end %i

%Scatter Solution (Q operator) to make it Local
Lmatrix=zeros(npoin,npoin);
Rvector=zeros(npoin,1);
icounter_CG=zeros(npoin_CG,1);
for i=1:npoin
    ip_CG=DG_to_CG(i); 
    icounter_CG(ip_CG)=icounter_CG(ip_CG) + 1;
    if icounter_CG(ip_CG) == 1
        Rvector(i)=Rvector_continuous(ip_CG);
        for j=1:npoin
            jp_CG=DG_to_CG(j);
            Lmatrix(i,j)=Lmatrix_continuous(ip_CG,jp_CG);
        end
    elseif icounter_CG(ip_CG) > 1
        Lmatrix(i,:)=0;
        Lmatrix(i,i)=1;
        Rvector(i)=0;
    end
end %i