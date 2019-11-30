%---------------------------------------------------------------------%
%This function updates/stores the frames of the movie for a CG/DG code.
%Written by F.X. Giraldo on September 16, 2014
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%Output Paramters: qi_movie(201,201,nframes+1) is the solution interpolated
%                  for plotting the movie
%                  time_movie(nframes+1,1) is the time for each frame
%                  xi and yi are the coordinates from INITIALIZE_MOVIE
%                  iframe is the new frame number
%Output Paramters: Movie_frames stores the actual movie in case you want to
%                  store it to a file.
%---------------------------------------------------------------------%
function Movie_frames = plot_movie(qi_movie,time_movie,xi,yi,iframe)

figure;
for i=1:iframe
    mesh(xi,yi,qi_movie(:,:,i));
    axis([-1 +1 -1 +1 0 1]);
    title_text=['Time = ' num2str(time_movie(i))];
    title([title_text],'FontSize',18);
    set(gca, 'FontSize', 18);
    colorbar('SouthOutside');
    M_i=getframe(gcf);
    Movie_frames(i)=M_i;
    pause(0.2);
end