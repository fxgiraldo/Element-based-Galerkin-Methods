%---------------------------------------------------------------------%
%A sample driver for how to use the Animation routines for Project 4
%MA4245
%Written by F.X. Giraldo on 9/17/2014
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

%To store or not to store...(the movie)
store_movie=0;
nframes=10; %Number of Frames in movie

%CREATE GRID, etc.

%Intialize Movie
[qi_movie,time_movie,xi,yi,iframe] = initialize_movie(coord,nframes);

%Time Integration
time=0;
for itime=1:ntime
   itime;
   time=time + dt;
  
   %BUILD CG/DG RHS, etc.
   
  %Evolve forward in Time
  for e=1:nelem
      qp(e,:,:)=a0*q0(e,:,:) + a1*q1(e,:,:) + dtt*rhs(e,:,:);
  end
   
   %Update/Store MOVIE frames
   [qi_movie,time_movie,iframe] = update_movie(qi_movie,time_movie,qp,coord,intma,npoin,nelem,ngl,ntime,itime,nframes,iframe,xi,yi,time);

end %itime

%Plot Movie
Movie_frames = plot_movie(qi_movie,time_movie,xi,yi,iframe);
if (store_movie == 1)
    movie2avi(M,'Movie.avi','fps',5);
end