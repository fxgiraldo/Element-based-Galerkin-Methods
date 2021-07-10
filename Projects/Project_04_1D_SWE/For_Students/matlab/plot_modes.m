%---------------------------------------------------------------------%
%This function plots the Modal Solution and Limiter
%Written by F.X. Giraldo on 1/2011
%           Naval Postgraduate School 
%           Monterey, CA 93943-5216
%---------------------------------------------------------------------%
function [test] = plot_modes(qmodal,Ipointer,limit_element,intma,coord,qp,nelem,ngl,nq,nop,noq,time,diss,main_text)

   test = 0;
   II=Ipointer';
   QQ(:,:)=qmodal';
   icheck=[1:ngl];

   %-------PLOT MODES------%
   subplot(3,1,1);
   %bar(II(1:nelem,icheck),QQ(1:nelem,icheck),'stacked');
   bar(QQ(1:nelem,icheck));
   %xlabel('Methods','FontSize',18);
   %ylabel('Modes','FontSize',18);
   %legend('0','1','2','3','4','5','6','7','8');
   title_text=[main_text ': diss = ', num2str(diss) ': Ne = ' num2str(nelem) ', N = ' num2str(nop) ', Q = ' num2str(noq) ', Time = ' num2str(time)];
   title([title_text],'FontSize',18);
   qmin=min(min(qmodal));
   qmax=max(max(qmodal));
   axis([ 1 nelem qmin qmax]);
   set(gca,'FontSize',18);
   %hold on;

   %-------PLOT Limiter------%
   subplot(3,1,2);
   bar(limit_element);
   qmin=min(limit_element);
   qmax=max(limit_element);
   qmin=0;
   qmax=1;
   axis([ 1 nelem qmin qmax]);
   %ylabel('Limiter','FontSize',18);
   set(gca,'FontSize',18);

   %-------PLOT Solution------%
   subplot(3,1,3)
   qq=zeros(nq,nelem);
    for ie=1:nelem
        for i=1:ngl
            ip=intma(i,ie);
            x(i)=coord(i,ie);
            y1(i)=qp(ip,1);
            y2(i)=qp(ip,2);
        end

        plot_handle=plot(x,y1,'r-',x,y2,'b:','LineWidth',2);
        %set(plot_handle,'LineWidth',2);
        hold on;
    end
%     axis([ -1 +1 -0.25 1.25]);
    %ylabel('Solution','FontSize',18);
    set(gca,'FontSize',18);
%        M_i=getframe(gcf);
%        M_modes(i)=M_i;
    pause(0.1);
   hold off;
   test=1;
