%-------------------------------------------------------
%This code applies a positivity limiter initially proposed
%by Xing, Zhang, and Shu in Positivity preserving high 
%order well balanced discontinuous Galerkin methods for 
%the shallow water equations (2010)
%-------------------------------------------------------
function [qp] = limiter_Shu_Positivity_Preserving(qp,nelem,psi,ngl,nq,wnq)
qlocal = zeros(ngl,ngl);
eps=1e-15;

for e=1:nelem
    
    %Calculate water height
    qlocal(:,:) = qp(e,:,:);
    
    %Loop Integration Points
    qe = zeros(nq,nq); %U for quad points
    for k2=1:nq
    for k1=1:nq
        q_k=0;
        %Interpolate at Integration Points
	    for i2=1:ngl
        for i1=1:ngl
            h_k=psi(i1,k1)*psi(i2,k2);
            q_k=q_k + h_k*qlocal(i1,i2);
        end %i1
        end %i2
        qe(k1,k2)=q_k;
    end %k1
    end %k2
    qmin  = min(qe(:));
    
    %Calculate Integral to find Mean Value
    qmean=0;
    for j=1:nq
    for i=1:nq
        qmean = qmean + wnq(i)*wnq(j)*qe(i,j);
    end
    end
    qmean=qmean/4;
    
    if (qmean < 0) 
        qmean
    end 
    
    %Calculate Theta
    theta = min(1, (qmean - eps)/( (qmean-qmin) ));

    for j=1:ngl
    for i=1:ngl
        qlocal(i,j) = (theta*(qlocal(i,j)-qmean)) + qmean;
        if (qmin == 0 & qmin == qmean)
            qlocal(i,j)= 0;
        end
    end
    end
    qp(e,:,:)=qlocal(:,:);
    
end


