%---------------------------------------------------------------------%
%This function plots the Modal Solution and Limiter
%Written by F.X. Giraldo on 1/2011
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [p_movie,u_movie,time_movie,mass_movie,energy_movie,L2_norm,iframe] = store_movie_frames(qp,q0,qb,intma,coord,wnq,psi,npoin,nelem,ngl,nq,p_movie,u_movie,time_movie,mass_movie,energy_movie,iframe,mass0,energy0,time,icase,eps)

iframe=iframe + 1;

p_movie(:,iframe)=qp(:,1);
u_movie(:,iframe)=qp(:,2);
time_movie(iframe)=time;

%Compute Mass
m1=0; m2=0;
for e=1:nelem
   dx=coord(ngl,e)-coord(1,e);
   jac=dx/2;
    for l=1:nq
      wq=wnq(l)*jac;
      for j=1:ngl
         jp=intma(j,e);
         h=qp(jp,1)+qb(jp);
         U=qp(jp,2);
         m1=m1 + wq*h*psi(j,l);
         m2=m2 + wq*(U)*psi(j,l);
      end
    end
end
mass_movie(iframe)=abs(m1-mass0);%/mass0;
energy_movie(iframe)=abs(m2-energy0);

%Compute Norm
[qe,qb,gravity] = exact_solution(intma,coord,npoin,nelem,ngl,time,icase,eps);
h_top=norm(q0(:,1)-qe(:,1),2);
h_bot=norm(qe(:,1),2);
u_top=norm(q0(:,2)-qe(:,2),2);
u_bot=norm(qe(:,2),2);
L2_norm(1,iframe)=h_top;   
L2_norm(2,iframe)=u_top;   


