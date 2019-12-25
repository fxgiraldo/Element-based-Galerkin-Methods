%---------------------------------------------------------------------%
%Written by M.A. Kopera
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function plot_grid_sub(coord,intma,sfc,nsfc,ngl)

%Plot Grid

    x=zeros(5,1);
    y=zeros(5,1);

    hold on;
   
    for e=1:nsfc
        el=sfc(e);
            for j=1:ngl-1
                for i=1:ngl-1
                    i1=intma(el,i,j);
                    i2=intma(el,i+1,j);
                    i3=intma(el,i+1,j+1);
                    i4=intma(el,i,j+1);
                    x(1)=coord(i1,1); y(1)=coord(i1,2);
                    x(2)=coord(i2,1); y(2)=coord(i2,2);
                    x(3)=coord(i3,1); y(3)=coord(i3,2);
                    x(4)=coord(i4,1); y(4)=coord(i4,2);
                    x(5)=coord(i1,1); y(5)=coord(i1,2);
                    plot_handle=plot(x,y,'-r');
                    set(plot_handle,'LineWidth',1.5);
                end
            end
             i1=intma(el,1,1);
             i2=intma(el,ngl,1);
             i3=intma(el,ngl,ngl);
             i4=intma(el,1,ngl);
             x(1)=coord(i1,1); y(1)=coord(i1,2);
             x(2)=coord(i2,1); y(2)=coord(i2,2);
             x(3)=coord(i3,1); y(3)=coord(i3,2);
             x(4)=coord(i4,1); y(4)=coord(i4,2);
             x(5)=coord(i1,1); y(5)=coord(i1,2);
             plot_handle=plot(x,y,'-b');
             set(plot_handle,'LineWidth',2);
    end

    %title_text=['Grid Plot For: Ne = ' num2str(nsfc) ', N = ' num2str(nop)];
%     title([title_text],'FontSize',18);      

    xlabel('X','FontSize',18);
    ylabel('Y','FontSize',18);
    set(gca,'FontSize',16);
    axis image
    



end
