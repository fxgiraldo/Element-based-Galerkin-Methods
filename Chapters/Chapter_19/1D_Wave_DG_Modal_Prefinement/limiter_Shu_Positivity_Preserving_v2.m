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
function [qpmod] = limiter_Shu_Positivity_Preserving_v2(qp,nelem,psi,nq,wnq,element_order)
qpmod = zeros(nq,nelem);
eps=1e-15;

for e=1:nelem
    N=element_order(e)+1;
    
    %Calculate water height
    qpmod(:,e)=0;
    qpmod(1:N,e) = qp(1:N,e);
    
    if (N>2)
    
        %--------------evaluate mj--------------------------
        Uq = zeros(nq,1); %U for quad points
        Uq = psi'*qpmod(:,e);
        qmin  = min(Uq(:));
        %--------------evaluate mj--------------------------

        %Calculate Integral for Ub
        qmean=0;
        for i=1:nq
            qmean = qmean + Uq(i)*wnq(i);
        end
        qmean=qmean/2; %size of computational element

        %Calculate Theta
        theta = min(1, (qmean - eps)/((qmean-qmin)));

        for i=2:N
            qpmod(i,e) = theta*(qpmod(i,e)-qmean) + qmean;
            if (qmin == 0 && qmin == qmean)
                qpmod(i,e)= 0;
            end
        end
        
    end %N>2
    
end

