%---------------------------------------------------------------------%
%This function plots the movie
%Written by F.X. Giraldo on 1/2011
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [test] = plot_solution(q0,qe,intma,coord,nelem,ngl,nop,noq,diss,time,l2_norm,store_plot,main_text)
test=0;
figure;

%Compute a gridpoint solution
np=nelem*ngl;
q_sol=zeros(np,1);
qe_sol=zeros(np,1);
x_sol=zeros(np,1);
ip=0;
for ie=1:nelem
for i=1:ngl
      ip=ip+1;
      q_sol(ip)=q_sol(ip) + q0(i,ie);
      qe_sol(ip)=qe_sol(ip) + qe(i,ie);
      x_sol(ip)=coord(i,ie);
end 
end

h=figure;
figure(h);
plot_handle=plot(x_sol,q_sol,'r-');
set(plot_handle,'LineWidth',2);
hold on
plot_handle=plot(x_sol,qe_sol,'b--');
set(plot_handle,'LineWidth',2);
axis([-1 +1 -0.25 +1.25]);

%Plot Discontinuous Elements
for ie=1:nelem
    x1=coord(1,ie);
    x2=coord(ngl,ie);
    y1=q0(1,ie);
    y2=q0(ngl,ie);
    plot(x1,y1,'ko');
    plot(x2,y2,'ko');
end

%Plot Continuous Elements
for ie=1:nelem
    i1=intma(1,ie);
    i2=intma(ngl,ie);
    x1=x_sol(i1);
    x2=x_sol(i2);
    y1=q_sol(i1);
    y2=q_sol(i2);
    plot(x1,y1,'ko');
    plot(x2,y2,'ko');
end

xlabel('x','FontSize',18);
ylabel('q(x,t)','FontSize',18);

% if (diss == 0)
%    file_ps=['dg_lgl_n' num2str(nelem) 'p' num2str(nop)];
%     legend('DG LGL','Exact');	
% elseif (diss == 1)
%    file_ps=['dg_lgl_upwind_n' num2str(nelem) 'p' num2str(nop)];
%    legend('DG LGL','Exact');	
% end

title_text=[main_text ': Diss = ' num2str(diss) ': Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', L2 Norm = ' num2str(l2_norm) ', T = ' num2str(time)];
title([title_text],'FontSize',18);
set(gca, 'FontSize', 18);

if store_plot == 1
    eval(['print ' file_ps ' -depsc']);
end
test=1;