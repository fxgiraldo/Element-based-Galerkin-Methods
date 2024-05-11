%---------------------------------------------------------------------%
%This function plots the Modal Solution and Limiter
%Written by F.X. Giraldo on 1/2011
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [p_movie,time_movie,mass_movie,L2_norm,iframe] = store_movie_frames(qp,q0,intma,coord,wnq,npoin,nelem,ngl,p_movie,time_movie,mass_movie,iframe,mass0,time,icase)

iframe=iframe + 1;

p_movie(:,iframe)=qp(:);
time_movie(iframe)=time;

[m1] = compute_Mass_and_Energy(qp,intma,coord,wnq,nelem,ngl);
mass_movie(iframe)=m1/mass0;

%Compute Norm
[qe] = exact_solution(intma,coord,npoin,nelem,ngl,time,icase);
h_top=norm(q0(:)-qe(:),2);
h_bot=norm(qe(:),2);
L2_norm(iframe)=h_top;   

