%---------------------------------------------------------------------%
%This function computes the RHS Flux contribution for the 1D SWE.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function rhs = create_rhs_flux(rhs,qp,qb,intma,nelem,ngl,diss,gravity,delta_nl,lsponge,nrbc_right,form_method)

%Initialize flag for Weak or Strong method
delta_form=0;
if strcmp(form_method,'strong')
    delta_form=1.0;
end

%Initialize
q_l=zeros(3,1);
q_r=zeros(3,1);

%-----------------------Left BC
er=1;
%Right Element
ir=intma(1,er);
hs_r=qp(ir,1);
U_r =qp(ir,2);
hb_r=qb(ir);
h_r =hs_r + hb_r;
u_r =U_r/(h_r);
%Left Element
h_l =h_r;
hs_l=hs_r;
U_l =-U_r; %hard wall boundaries
hb_l=hb_r;
u_l =U_l/(h_l);
%Store States
q_l(1)=hs_l;
q_l(2)=hb_l;
q_l(3)=u_l;
q_r(1)=hs_r;
q_r(2)=hb_r;
q_r(3)=u_r;
%Rusanov Flux
[flux,flux_l,flux_r] = rusanov_flux(q_l,q_r,diss,gravity,delta_nl);
rhs(ir,1)=rhs(ir,1) + (flux(1)-delta_form*flux_r(1));
rhs(ir,2)=rhs(ir,2) + (flux(2)-delta_form*flux_r(2));

%-----------------------Right BC
%if (lsponge==0)
    el=nelem;
    %Left Element
    il=intma(ngl(el),el);
    hs_l=qp(il,1);
    U_l =qp(il,2);
    hb_l=qb(il);
    h_l =hs_l + hb_l;
    u_l =U_l/(h_l);
    %Right Element
    h_r =h_l;
    hs_r=hs_l;
    if (nrbc_right==1)
        U_r =-U_l; %hard wall boundary=>NFBC
    elseif (nrbc_right==2)
        U_r =U_l; %extrapolation/outflow boundary=>Neumann
    end
    hb_r=hb_l;
    u_r =U_r/(h_r);   
    %Store States
    q_l(1)=hs_l;
    q_l(2)=hb_l;
    q_l(3)=u_l;
    q_r(1)=hs_r;
    q_r(2)=hb_r;
    q_r(3)=u_r;
    %Rusanov Flux
    [flux,flux_l,flux_r] = rusanov_flux(q_l,q_r,diss,gravity,delta_nl);
    rhs(il,1)=rhs(il,1) - (flux(1)-delta_form*flux_l(1));
    rhs(il,2)=rhs(il,2) - (flux(2)-delta_form*flux_l(2));
%end