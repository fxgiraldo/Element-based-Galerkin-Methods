%----------------------------------------------------------------------%
%This subroutine builds the FLUX vector for the Strong Form DGM-SEM
%on Quadrilateral Elements for the 2D Euler Equations.
%Written by Francis X. Giraldo on 1/2001
%           Naval Postgraduate School
%           Department of Applied Mathematics
%           Monterey, CA 93943-5216
%----------------------------------------------------------------------%
function rhs = compute_flux_dg_TensorProduct_inexact_1(rhs,q,u,v,nx,ny,jac_side,...
    ngl,imapl,imapr,ffc,nffc,face,P1g,P2g,P1s,P2s)

%local arrays
ql=zeros(ngl,1);
qr=zeros(ngl,1);
ul=zeros(ngl,1);
vl=zeros(ngl,1);

flux_qla=zeros(ngl,1);
flux_qr1=zeros(ngl,1);
flux_qr2=zeros(ngl,1);
flux_qp=zeros(ngl,1);
flux_qc1=zeros(ngl,1);
flux_qc2=zeros(ngl,1);

%Construct FVM-type Operators
for ifa=1:nffc
    is=ffc(ifa);

    if face(is,9)==0
       %Store Left Side Variables
       iel=face(is,5);%psideh(is,3);

       if (face(is,10) ~= -6) %periodic bc
          ilocl=face(is,3);
          for l=1:ngl
              %Get Pointers
              il=imapl(ilocl,1,l);
              jl=imapl(ilocl,2,l);
              %Left Element
              ql(l)=q(iel,il,jl);
              ul(l)=u(iel,il,jl);
              vl(l)=v(iel,il,jl);
          end %l  

          %Store Right Side Variables
          ier=face(is,6);
          if (ier > 0 ) 
             ilocr=face(is,4);            
             for l=1:ngl
                 %Get Pointers
                 ir=imapr(ilocr,1,l);
                 jr=imapr(ilocr,2,l);
                 %Right Element
                 qr(l)=q(ier,ir,jr);
             end              %l  
          end %if ier
            
          flux_q = rusanov_flux(jac_side,nx,ny,ql,qr,ul,vl,ngl,is);
          
          for l=1:ngl
              wq=jac_side(is,l);
              %--------------Left Side------------------%
              il=imapl(ilocl,1,l);
              jl=imapl(ilocl,2,l);
              rhs(iel,il,jl)=rhs(iel,il,jl) - wq*flux_q(l);

              %--------------Right Side------------------%
              if (ier > 0)
                 ir=imapr(ilocr,1,l);
                 jr=imapr(ilocr,2,l);
                 rhs(ier,ir,jr)=rhs(ier,ir,jr) + wq*flux_q(l);
              end %if 

          end %l   
          
       end %if iel face  
    elseif face(is,9)==1

        %find parent and children elements
         if face(is,7)==0
             parent = face(is,5);
             child1 = face(is,6);
             child2 = face(is,8);
             ipf = 3;
             icf = 4;
%              sg = -1;
         elseif face(is,8)==0
             parent = face(is,6);
             child1 = face(is,5);
             child2 = face(is,7);
             ipf = 4;
             icf = 3;
%              sg=1;
         end
             
        
       %Store Parent Variables

       if (face(is,10) ~= -6) %periodic bc
          ilocl=face(is,ipf);
          for l=1:ngl
              %Get Pointers
              il=imapl(ilocl,1,l);
              jl=imapl(ilocl,2,l);
              %Parent Element
              ql(l)=q(parent,il,jl);
              ul(l)=u(parent,il,jl);
              vl(l)=v(parent,il,jl);
          end %l  

          %Store Children Variables

             ilocr=face(is,icf);            
             for l=1:ngl
                 %Get Pointers
                 ir=imapr(ilocr,1,l);
                 jr=imapr(ilocr,2,l);
                 %Children Element
                 qc1(l)=q(child1,ir,jr);
                 qc2(l)=q(child2,ir,jr);
                 uc1(l)=u(child1,ir,jr);
                 uc2(l)=u(child2,ir,jr);
                 vc1(l)=v(child1,ir,jr);
                 vc2(l)=v(child2,ir,jr);
                 
             end              %l               
             
          %gather the data
%           [qla]=gather(qc1,qc2,ngl);  
          
          %gather scatter
%            [qla,qc1a,qc2a] = scatter_gather_l2(ql,qc1,qc2,ngl,P1,P2);
%          
%           qla = qc1*P1g+qc2*P2g;
          qc1a = ql'*P1s;
          qc2a = ql'*P2s;

          % compute Rusanov fluxes
%           flux_q = rusanov_flux(jac_side,nx,ny,ql,qla,ul,vl,ngl,is);
          flux_q1 = rusanov_flux(jac_side,nx,ny,qc1a,qc1,uc1,vc1,ngl,is);
          flux_q2 = rusanov_flux(jac_side,nx,ny,qc2a,qc2,uc2,vc2,ngl,is);
          
%           [flux_q1,flux_q2]=scatter(flux_q,ngl);
          
          flux_q = flux_q1'*P1g+flux_q2'*P2g;

          %apply flux
          for l=1:ngl
              wq=jac_side(is,l);
              wqc=wq/2;            
              
              %--------------Left Side------------------%
              il=imapl(ilocl,1,l);
              jl=imapl(ilocl,2,l);
%               rhs(parent,il,jl)=rhs(parent,il,jl) + sg*wq*flux_q(l);
%  
              rhs(parent,il,jl)=rhs(parent,il,jl) - wq*flux_q(l);
 

              %--------------Right Side------------------%
              ir=imapr(ilocr,1,l);
              jr=imapr(ilocr,2,l);
%               rhs(child1,ir,jr)=rhs(child1,ir,jr) - sg*wqc*flux_q1(l); %sg is a +/- sign depending on face set-up
%               rhs(child2,ir,jr)=rhs(child2,ir,jr) - sg*wqc*flux_q2(l);  
               rhs(child1,ir,jr)=rhs(child1,ir,jr) + wqc*flux_q1(l); %sg is a +/- sign depending on face set-up
              rhs(child2,ir,jr)=rhs(child2,ir,jr) + wqc*flux_q2(l);  
           end
           
       end
        
    end
end %is
end

