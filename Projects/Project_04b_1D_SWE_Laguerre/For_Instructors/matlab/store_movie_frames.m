%---------------------------------------------------------------------%
%This function plots the Modal Solution and Limiter
%Written by F.X. Giraldo on 1/2011
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [p_movie,u_movie,time_movie,mass_movie,energy_movie,iframe] = store_movie_frames(qp,qb,intma,jac,wnq,nelem,ngl,p_movie,u_movie,time_movie,mass_movie,energy_movie,iframe,mass0,energy0,time)

iframe=iframe + 1;

p_movie(:,iframe)=qp(:,1);
u_movie(:,iframe)=qp(:,2);
time_movie(iframe)=time;

%Compute Mass
[mass,energy]=compute_mass(qp,qb,intma,ngl,nelem,wnq,jac);
mass_movie(iframe)=abs(mass-mass0);%/mass0;
energy_movie(iframe)=abs(energy-energy0);