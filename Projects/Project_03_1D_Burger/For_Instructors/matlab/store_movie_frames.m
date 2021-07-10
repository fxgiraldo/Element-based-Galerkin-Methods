%---------------------------------------------------------------------%
%This function plots the Modal Solution and Limiter
%Written by F.X. Giraldo on 1/2011
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [p_movie,u_movie,time_movie,mass_movie,energy_movie,L2_norm,iframe] = store_movie_frames(qp,q0,intma,coord,wnq,psi,npoin,nelem,ngl,nq,p_movie,u_movie,time_movie,mass_movie,energy_movie,iframe,mass0,energy0,time,icase)

iframe=iframe + 1;

p_movie(:,iframe)=qp(:,1);
u_movie(:,iframe)=qp(:,2);
time_movie(iframe)=time;

[m1,m2] = compute_Mass_and_Energy(qp,intma,coord,wnq,psi,nelem,ngl,nq);
mass_movie(iframe)=m1/mass0;
energy_movie(iframe)=m2/energy0;

%Compute Norm
[qe] = exact_solution_dg(intma,coord,npoin,nelem,ngl,time,icase);
h_top=norm(q0(:,1)-qe(:,1),2);
h_bot=norm(qe(:,1),2);
u_top=norm(q0(:,2)-qe(:,2),2);
u_bot=norm(qe(:,2),2);
L2_norm(1,iframe)=h_top;   
L2_norm(2,iframe)=u_top;   


