%---------------------------------------------------------------------%
%This function initializes the movie frames for a CG/DG code.
%Written by F.X. Giraldo on September 16, 2014
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%Input Parameters: coord(npoin,2) are the coordinates from CREATE_GRID_2D
%                  nframes are the number of frames you want to store
%                  (user's choice).
%Output Paramters: qi_movie(201,201,nframes+1) is the solution interpolated
%                  for plotting the movie
%                  time_movie(nframes+1,1) is the time for each frame
%                  xi and yi are the coordinates reshaped for the movie
%                  iframe is initialized to iframe=0
%---------------------------------------------------------------------%
function [qi_movie,time_movie,xi,yi,iframe] = initialize_movie(coord,nframes)

iframe=0;
xmin=min(coord(:,1)); xmax=max(coord(:,1));
ymin=min(coord(:,2)); ymax=max(coord(:,2));
nxx=200; nyy=200;
qi_movie=zeros(nxx+1,nyy+1,nframes+1);
time_movie=zeros(nframes+1,1);
dx=(xmax-xmin)/nxx;
dy=(ymax-ymin)/nyy;
[xi,yi]=meshgrid(xmin:dx:xmax,ymin:dy:ymax);
    