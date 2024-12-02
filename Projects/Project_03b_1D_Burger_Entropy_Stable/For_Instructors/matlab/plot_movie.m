%---------------------------------------------------------------------%
%This function plots the movie
%Written by F.X. Giraldo on 1/2011
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [test] = plot_movie(qi_movie,time_movie,iframe,coord,nelem,ngl,nop,noq,diss,store_movie,main_text)

test=0;
figure;
for i=1:iframe

    for e=1:nelem
        for j=1:ngl
            x(j)=coord(j,e);
            y(j)=qi_movie(j,e,i);
        end
        plot_handle=plot(x,y,'r-');
        set(plot_handle,'LineWidth',2);
        hold on;
    end
    axis([-1 +1 -0.25 +1.25]);
    title_text=[main_text ': Diss = ' num2str(diss) ', Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time_movie(i))];
    title([title_text],'FontSize',18);
    hold on;

    %Plot Elements
    for e=1:nelem
        x1=coord(1,e);
        x2=coord(ngl,e);
        y1=qi_movie(1,e,i);
        y2=qi_movie(ngl,e,i);
        plot(x1,y1,'ko');
        plot(x2,y2,'ko');
    end

    title([title_text],'FontSize',18);
    xlabel('x','FontSize',18);
    ylabel('q(x,t)','FontSize',18);
    set(gca, 'FontSize', 18);
    M_i=getframe(gcf);
    M(i)=M_i;
    pause(0.1);
    hold off;
end
if store_movie == 1
    file_movie=['DG_LGL_n' num2str(nelem) '_p' num2str(nop) '_q' num2str(noq) '_diss' num2str(diss) '.avi'];
    movie2avi(M,file_movie,'fps',5);
end  
test=1;