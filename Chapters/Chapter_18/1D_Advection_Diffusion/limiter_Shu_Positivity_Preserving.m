%-------------------------------------------------------
%This code applies a positivity limiter initially proposed
%by Xing, Zhang, and Shu in ?Positivity preserving high 
%order well balanced discontinuous Galerkin methods for 
%the shallow water equations? (2010)
%
%Written by Haley Lane (NREIP) 6/2013-8/2013
%           Department of Applied Mathemcatics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%-------------------------------------------------------
function [qpmod] = limiter_Shu_Positivity_Preserving(qp, nelem, psi,ngl,nq,wnq)
qpmod = zeros(ngl,nelem);
eps=1e-15;

for e=1:nelem
    Uq = zeros(nq,1); %U for quad points
    qmean = 0;  %average U
    
    %Calculate water height
    qpmod(:,e) = qp(:,e);
    
    %--------------evaluate mj--------------------------
    Uq = psi'*qpmod(:,e);
    qmin  = min(Uq(:));
    %--------------evaluate mj--------------------------
    
    %Calculate Integral for Ub
    for i=1:nq
        qmean = qmean + Uq(i)*wnq(i);
    end
    qmean=qmean/2;
    
    %Calculate Theta
%     theta = min(1, Ub(1)/(abs(Ub(1)-m)+epsil));
    theta = min(1, (qmean - eps)/((qmean-qmin)));

    for i=1:ngl
        qpmod(i,e) = (theta*(qpmod(i,e)-qmean))+qmean;
        if (qmin == 0 & qmin == qmean)
            qpmod(i,e)= 0;
        end
    end

end

