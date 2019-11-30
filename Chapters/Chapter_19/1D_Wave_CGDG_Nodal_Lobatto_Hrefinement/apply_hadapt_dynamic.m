%---------------------------------------------------------------------%
%This function creates the H-refinement arrays for all levels of
%refinement.
%Written by F.X. Giraldo on 3/10/2016
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [qp,active] = apply_hadapt_dynamic(qp,active,ref_level,children,parent,hadapt_lev,sfc,nsfc,VDM_Inv_Matrix,nelem,ngl,refine_tol,coarsen_tol,P1g,P2g,P1s,P2s)

eps=1e-15;
qm=zeros(ngl,1);
q_parent=zeros(ngl,1);
q_child1=zeros(ngl,1);
q_child2=zeros(ngl,1);
N=ngl;

coarsen=zeros(nelem,1);
refine=zeros(nelem,1);

%Mark elements for Refinement/derefinement 
for ee=1:nsfc
    e=sfc(ee);

    %Get Modal Coefficients
    qm(1:N)=0;
    i=N;
    qm(1:i)=VDM_Inv_Matrix(1:i,1:i,i-1)*qp(1:i,e);

    %Check Criterion
    numerator=norm(qm(i:N));
    denominator=norm(qm(1:N))+eps;
    check=numerator/denominator;

    %Check if current mode is significant
    if (check <= coarsen_tol || denominator <= coarsen_tol) %Coarsen
       if ( ref_level(e) >= 1 )
           iparent=parent(e);
           coarsen(iparent)=coarsen(iparent) + 1;
       end %if
    elseif (check > refine_tol && denominator > refine_tol) %Refine
       if ( ref_level(e) < hadapt_lev )
           refine(e)=1;
       end %if
    end %if      
end %ee

%Apply Derefinement 
for e=1:nelem
    if (coarsen(e) == 2) 
        iparent=e;
        ichild1=children(1,e);
        ichild2=children(2,e);
        q_child1(:)=qp(:,ichild1);
        q_child2(:)=qp(:,ichild2);
        q_parent(:)=P1g*q_child1 + P2g*q_child2;
        
        %Modify Data
        qp(:,iparent)=q_parent(:);
        active(iparent)=1;
        active(ichild1)=0;
        active(ichild2)=0;
    end   
end %e

%Apply Refinement 
for ee=1:nsfc
    e=sfc(ee);
    if (refine(e) == 1) 
        iparent=e;
        ichild1=children(1,e);
        ichild2=children(2,e);
        q_parent(:)=qp(:,e);       
        q_child1(:)=P1s*q_parent(:);
        q_child2(:)=P2s*q_parent(:);
     
        %Modify Data
        qp(:,ichild1)=q_child1(:);
        qp(:,ichild2)=q_child2(:);
        active(iparent)=0;
        active(ichild1)=1;
        active(ichild2)=1;
    end   
end %e