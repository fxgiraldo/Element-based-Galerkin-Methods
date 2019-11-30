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
function qp = limiter_Shu_Positivity_Preserving(qp,nelem,element_order,nop_max,wgl_matrix,Interpolation_Matrix)
ngl=nop_max+1;
q=zeros(ngl,1);
eps=1e-15;

for e=1:nelem
    N=element_order(e)+1;
       
    %Get Nodal Coefficients
    q(:)=Interpolation_Matrix(:,1:N,N-1)*qp(1:N,e);

    %Min Value
    qmin  = min(qp(1:N,e));

    %Mean Value
    qmean=0;
    for k=1:ngl
        qmean = qmean + q(k)*wgl_matrix(k,nop_max);
    end
    qmean=qmean/2; %size of computational element

    %Calculate Theta
    theta = min(1, (qmean)/((qmean-qmin)));

    for k=1:N
        qp(k,e) = theta*(qp(k,e)-qmean) + qmean;
        if (qmin == 0 && qmean == 0)
           qp(k,e)= 0;
        end
        qp(k,e)=max(0,qp(k,e));
    end
end

