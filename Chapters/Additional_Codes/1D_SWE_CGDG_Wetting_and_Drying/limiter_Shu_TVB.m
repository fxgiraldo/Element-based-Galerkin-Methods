%-------------------------------------------------------
%This code applies a TVB limiter proposed
%by Shu in "Positivity preserving high 
%order well balanced discontinuous Galerkin methods for 
%the shallow water equations" (2010)
%
%Written by Shiva Gopalakrishnan 
%Modified by F.X. Giraldo to go into CGDG Unified code.
%           Department of Applied Mathemcatics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%-------------------------------------------------------
function [qp] = limiter_Shu_TVB(qp,qb,coord,vdm,nelem,ngl)
qp_avg=zeros(1,nelem);        
q_total=zeros(ngl,nelem);
m=1e-5;

for e=1:nelem
    for i=1:ngl
        q_total(i,e)=qp(1,i,e);% + qb(i,e);
    end
end

for e=1:nelem
    qp_modal=zeros(ngl);
    for i=1:ngl
        for j=1:ngl
            qp_modal(i) = qp_modal(i) + vdm(i,j)*q_total(j,e);
        end
    end
    qp_avg(1,e)=qp_modal(1);
end

for e=2:nelem-1
    qp_modal=zeros(ngl);
    iel=e-1;
    ier=e+1;
    if(e == 1)
        iel=nelem;
    end
    if(e == nelem)
        ier=1;
    end
    qp_1=q_total(1,e) - qp_avg(1,e);
    qp_2=-(q_total(ngl,e) - qp_avg(1,e));
    delta_plus=qp_avg(1,ier)-qp_avg(1,e);
    delta_minus=qp_avg(1,e)-qp_avg(1,iel);
    dx=coord(ngl,e) - coord(1,e);
    qp_1mod=minmod(qp_1,delta_plus,delta_minus,dx,m);
    qp_2mod=minmod(qp_2,delta_plus,delta_minus,dx,m);
    d1=qp_1 - qp_1mod;
    d2=qp_2 - qp_2mod;
    qp_nodal=zeros(ngl);
    qp_nodal(1)=qp_avg(1,e)+qp_1mod;
    qp_nodal(ngl)=qp_avg(1,e)-qp_2mod;
    for i=1:ngl
        qp(1,i,e)=qp_nodal(i);%-qb(i,e);
    end
end


