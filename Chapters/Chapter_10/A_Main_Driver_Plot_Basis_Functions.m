%---------------------------------------------------------------------%
%This code plots the basis functions.
%Written by F.X. Giraldo on 10/2003
%           Department of Applied Mathematics
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
clear all; 
close all;

%Input Data
nop=2;    %Interpolation Order
Ns=51;
plot_basis=1;
plot_grid=1;
store_basis=1;
store_grid=1;
point_method=1; %=1 LGL, =2 LG
ngl=nop + 1;

%Compute LGL Points
if point_method == 1
    [xgl,wgl]=legendre_gauss_lobatto(ngl); %LGL Points
elseif point_method == 2
    [xgl,wgl] = legendre_gauss(ngl);       %LG Points
end

%Compute Legendre Cardinal functions and derivatives
xs=linspace(-1,+1,Ns);
[psi,dpsi] = lagrange_basis2(ngl,xgl,Ns,xs);
[psi2,dpsi2] = lagrange_basis2(ngl,xgl,ngl,xgl);

if plot_basis == 1
    for j=1:ngl
        for i=1:ngl
            I=(j-1)*ngl + i
            for l=1:Ns
                for k=1:Ns
                    psis(k,l)=psi(i,k)*psi(j,l);
                end
            end
            psis=psis';

            %Plot Basis Function Values
            figure;
            surf(xs,xs,psis);
            xlabel('\xi','FontSize',18);
            ylabel('\eta','FontSize',18);
            zlabel('\psi(\xi,\eta)','FontSize',18);
            axis image
            title_text=[' N = ' num2str(nop),' \psi', num2str(I)];
            if store_basis == 0
                title([title_text],'FontSize',18);
            end
            set(gca, 'FontSize', 18);
            axis([ -1 1 -1 1 -1 2.5]);
            axis image
            %whitebg('black');
            hold on;

            %Plot Nodes
            ip=0;
            for l=1:ngl
                for k=1:ngl
                    ip=ip+1;
                    xx(ip)=xgl(k);
                    yy(ip)=xgl(l);
                    psis2(ip)=psi2(i,k)*psi2(j,l);
                end
            end
            plot3(xx,yy,psis2,'ro','LineWidth',2);
            if point_method == 1
                file_ps=['LGL_Quadrilateral_Grid_N' num2str(nop) '_psi_' num2str(I)];
            elseif point_method == 2
                file_ps=['LG_Quadrilateral_Grid_N' num2str(nop) '_psi_' num2str(I)];
            end
            if store_basis == 1
                eval(['print ' file_ps ' -depsc']);
            end

        end%i
    end %j
end

if (plot_grid == 1)
    x=zeros(5,1);
    y=zeros(5,1);
    figure;
    hold on;
    
    for j=1:ngl-1
        j1=j;
        j2=j+1;
        for i=1:ngl-1
            i1=i;
            i2=i+1;
            x(1)=xgl(i1); y(1)=xgl(j1);
            x(2)=xgl(i2); y(2)=xgl(j1);
            x(3)=xgl(i2); y(3)=xgl(j2);
            x(4)=xgl(i1); y(4)=xgl(j2);
            x(5)=xgl(i1); y(5)=xgl(j1);
            plot_handle=plot(x,y,'-ro');
            set(plot_handle,'LineWidth',2);
        end
    end
    title_text=['Grid Plot For:  N = ' num2str(nop)];
    if store_grid == 0
        title([title_text],'FontSize',18);      
    end

    xlabel('\xi','FontSize',18);
    ylabel('\eta','FontSize',18);
    set(gca,'xtick',[-1:0.5:1]);
    set(gca,'ytick',[-1:0.5:1]);
    set(gca, 'FontSize', 18);
    axis([ -1 1 -1 1]);
    axis square
    axis image
    
    if point_method == 1
        file_ps=['LGL_Quadrilateral_Grid_N' num2str(nop)];
    elseif point_method == 2
        file_ps=['LG_Quadrilateral_Grid_N' num2str(nop)];
    end
    if store_grid == 1
       eval(['print ' file_ps ' -depsc']);
    end

end
    