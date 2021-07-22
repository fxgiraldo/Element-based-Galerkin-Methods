%---------------------------------------------------------------------%
%This function computes the Derivative Mapping for General 2D Quad Grids.
%Written by F.X. Giraldo on 4/2008
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [f_ksi,f_eta,f_zeta] = map_deriv(psi,dpsi,f,ngl,nq)

%initialize
f_ksi=zeros(nq,nq,nq);
f_eta=zeros(nq,nq,nq);
f_zeta=zeros(nq,nq,nq);


for n=1:nq
for m=1:nq
for l=1:nq
   
    sum_ksi=0;
    sum_eta=0;
    sum_zeta=0;

    for k=1:ngl
    for j=1:ngl
    for i=1:ngl
        sum_ksi= sum_ksi  + dpsi(i,l)*psi(j,m)*psi(k,n)*f(i,j,k);
        sum_eta =sum_eta  + psi(i,l)*dpsi(j,m)*psi(k,n)*f(i,j,k);
        sum_zeta=sum_zeta + psi(i,l)*psi(j,m)*dpsi(k,n)*f(i,j,k);    
    end %i
    end %j
    end %k
   
    f_ksi(l,m,n)=sum_ksi;
    f_eta(l,m,n)=sum_eta;
    f_zeta(l,m,n)=sum_zeta;


end %l
end %m
end %n
