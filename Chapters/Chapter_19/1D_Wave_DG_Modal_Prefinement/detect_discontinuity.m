%---------------------------------------------------------------------%
%This function computes the LGL grid and elements.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [limit_element,qp] = detect_discontinuity(qp,nelem,ngl,nq,L,Ls,dx,xmu)

%Initialize
limit_element=zeros(nelem,1);
dx_p=dx^(ngl/2);
tol=1;

%Integrate Flux Terms
for ie=1:nelem
   iel=ie-1;
   ier=ie;
   if (ie == 1) 
      iel=nelem;
   end 
   
   %Jacobians
   jac=dx/2;
   ksi_x=2/dx;
   
   %Interpolate F (via q) onto element side/interface/edges
   q_l=0; q_r=0;
   for i=1:ngl
       q_l=q_l + Ls(i,2)*qp(i,iel);
       q_r=q_r + Ls(i,1)*qp(i,ier);
   end
   
   %Compute Average Q at Integration Points
   q_average=0;
   for k=1:nq
      q_k=0;
      for j=1:ngl
          q_k=q_k + L(j,k)*qp(j,ie);
      end
      q_average=q_average + q_k/nq;
   end %l
   
   %Icheck=abs(q_r - q_l)/(q_average);
   %Icheck=abs(q_r - q_l);
   Icheck=abs(q_r - q_l)/(dx_p);
   %limit_element(ie)=Icheck;
   if Icheck > tol
       limit_element(ie)=1;
   end
   
   %Limit Solution
   if limit_element(ie) == 1 
    ibeg=1;
    iend=ngl;
    for i=1:ibeg
    weight(i)=1;
    end
    for i=ibeg+1:iend
        weight(i)=1.0 - ((i-ibeg)/(iend-ibeg))^2.5;
    end
    
    qp(:,ie)= weight(:).*qp(:,ie);
   end %if
   
end %ie
