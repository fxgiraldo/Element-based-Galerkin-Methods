%---------------------------------------------------------------------%
%This code plots the Stability of Time-Differencing Schemes.
%Written by F.X. Giraldo on 10/07
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%Variables:
%---------------------------------------------------------------------%

%clear; clear all; clf; close all;
clear all
close all

%Read Number of Files
iplot=0;
% delta=input(' 0=Explicit \n 1=Implicit \n Enter delta: ');
method=input(' 1=BDF1 \n 2=BDF2 \n 3=BDF3 \n 4=BDF4 \n 5=BDF5 \n 6=BDF6 \n 7=AM2/AI2 \n 8=LF2 \n Enter Method: ');
   
%setup TDS coefficients
if (method <= 6)
    text=['BDF' num2str(method)];
elseif (method == 7)
    text=['AI2'];
elseif (method == 8)
    text=['LF2'];
end

ns=500;
nf=500;
csmax=2;
cfmax=20;
ksdtmin=-csmax;
ksdtmax=+csmax;
dksdt=(ksdtmax-ksdtmin)/(ns-1);
kfdtmin=0;
kfdtmax=+cfmax;
dkfdt=(kfdtmax-kfdtmin)/(nf-1);
amp=zeros(2,ns,nf);
kfdt=zeros(nf,1);
ksdt=zeros(ns,1);

if (method == 1)
    %BDF1
    K=1;
    a(1)=1;
    gamma=1;
    b(1)=1;
    rp=1;
    r(1)=0;
elseif (method == 2)
    %BDF2
    K=2;
    a(1)=4/3;
    a(2)=-1/3;
    gamma=2/3;
    b(1)=2;
    b(2)=-1;
    rp=1;
    r(1)=0;
    r(2)=0;
elseif (method == 3)
   %BDF3
   K=3;
   a(1)=18/11;
   a(2)=-9/11;
   a(3)=2/11;
   gamma=6/11;
   b(1)=3;
   b(2)=-3;
   b(3)=1;
   rp=1;
   r(1)=0;
   r(2)=0;
   r(3)=0;
elseif (method == 4)
   %BDF4
   K=4;
   a(1)=48/25;
   a(2)=-36/25;
   a(3)=16/25;
   a(4)=-3/25;
   gamma=12/25;
   b(1)=4;
   b(2)=-6;
   b(3)=4;
   b(4)=-1;
   rp=1;
   r(1)=0;
   r(2)=0;
   r(3)=0;
   r(4)=0;
elseif (method == 5)
   %BDF5
   K=5;
   a(1)=300/137;
   a(2)=-300/137;
   a(3)=200/137;
   a(4)=-75/137;
   a(5)=12/137;
   gamma=60/137;
   b(1)=5;
   b(2)=-10;
   b(3)=10;
   b(4)=-5;
   b(5)=1;
   rp=1;
   r(1)=0;
   r(2)=0;
   r(3)=0;
   r(4)=0;
   r(5)=0;
elseif (method == 6)
   %BDF6
   K=6;
   a(1)=360/147;
   a(2)=-450/147;
   a(3)=400/147;
   a(4)=-225/147;
   a(5)=72/147;
   a(6)=-10/147;
   gamma=60/147;
   b(1)=6;
   b(2)=-15;
   b(3)=20;
   b(4)=-15;
   b(5)=6;
   b(6)=-1;
   rp=1;
   r(1)=0;
   r(2)=0;
   r(3)=0;
   r(4)=0;
   r(5)=0;
   r(6)=0;
elseif (method == 7)
   %AM2/AI2
   K=3;
   a(1)=1;
   a(2)=0;
   a(3)=0;
   gamma=1;
   b(1)=23.0/12.0;
   b(2)=-16.0/12.0;
   b(3)=5.0/12.0;
   rp=15.0/12.0;
%    r(1)=-35.0/12.0; %This is what is in NUMA because
%    r(2)=25.0/12.0;  %in Implicit part we need [ F-{S} ]
%    r(3)=-5.0/12.0;
   r(1)=-12.0/12.0;
   r(2)=9.0/12.0;
   r(3)=0;
elseif (method == 8)
   alpha=input(' Enter LF2 Implicit weighting (0,1): ');
   %LF2
   K=2;
   a(1)=0;
   a(2)=1;
   gamma=2;
   b(1)=1;
   b(2)=0;
   rp=alpha;
   r(1)=0;
   r(2)=(1.0-alpha);
end

c=zeros(K+1,1);
for l=1:nf
   kfdt(l)=kfdtmin + (l-1)*dkfdt;
   for k=1:ns
      ksdt(k)=ksdtmin + (k-1)*dksdt;
      c(1)=1 - gamma*1i*kfdt(l)*rp;
      for m=2:K+1
	      c(m)=-(a(m-1) + gamma*b(m-1)*1i*ksdt(k) + gamma*r(m-1)*1i*kfdt(l) );
      end
      rr=roots(c);
      for m=1:K
           amp(m,k,l)=sqrt( rr(m)*conj(rr(m)) );
      end
   end       
end


[X,Y]=meshgrid(ksdt,kfdt);
V=0:0.05:1.05;
amp1=zeros(ns,nf);
amp2=zeros(ns,nf);
amp_max=zeros(ns,nf);
amp_min=zeros(ns,nf);
for l=1:nf
    for k=1:ns
         amp_max(k,l)=max(amp(:,k,l));
         amp_min(k,l)=min(amp(:,k,l));
%          if (amp_max(k,l) > 1 + 1e-6)
         if (amp_max(k,l) > 1 + 1e-3)
             amp_max(k,l) = nan;
         end
         for m=1:K
%              if (amp(m,k,l) > 1 + 1e-6)
             if (amp(m,k,l) > 1 + 1e-3)  
                 amp(m,k,l) = nan;
             end
         end
    end
end


% h=figure;
% [cs,h]=contourf(X',Y',amp_max,V);
% hold on
% set(gca, 'CLim', [0,1]);
% colorbar('SouthOutside');
% set(gca,'FontSize',18);
% xlabel('k_s \Delta t','FontSize',18); 
% ylabel('k_f \Delta t','FontSize',18); 
% %set(gca,'XTick',[ksdtmin:1:ksdtmax]);
% %set(gca,'YTick',[kfdtmin:1:kfdtmax]);
% axis([ksdtmin ksdtmax kfdtmin kfdtmax]);
% %axis square
% %title(text);
% 
% h=figure;
% [cs,h]=contourf(X',Y',amp_min,V);
% hold on
% set(gca, 'CLim', [0,1]);
% colorbar('SouthOutside');
% set(gca,'FontSize',18);
% xlabel('k_s \Delta t','FontSize',18); 
% ylabel('k_f \Delta t','FontSize',18); 
% %set(gca,'XTick',[ksdtmin:1:ksdtmax]);
% %set(gca,'YTick',[kfdtmin:1:kfdtmax]);
% axis([ksdtmin ksdtmax kfdtmin kfdtmax]);
% %axis square
% %title(text);

for i=1:K
    h=figure;
    amp_i=squeeze(amp(i,:,:));
    [cs,h]=contourf(X',Y',amp_i,V);
    hold on
    set(gca, 'CLim', [0,1]);
    colorbar('SouthOutside');
    colormap(jet);
    set(gca,'FontSize',18);
    xlabel('k_s \Delta t','FontSize',18); 
    ylabel('k_f \Delta t','FontSize',18); 
    %set(gca,'XTick',[ksdtmin:1:ksdtmax]);
    %set(gca,'YTick',[kfdtmin:1:kfdtmax]);
    axis([ksdtmin ksdtmax kfdtmin kfdtmax]);
    %axis square
    title([ ' K = ',num2str(i)]);
end

% figure
% subplot(1,2,1);
% plot(ksdt,amp_max(:,1),'r-','LineWidth',2);
% hold on;
% plot(ksdt,amp_min(:,1),'b--','LineWidth',2);
% set(gca,'FontSize',18);
% xlabel('k_s \Delta t','FontSize',18); 
% ylabel('|A|','FontSize',18); 
% axis([ksdtmin ksdtmax 0 1.75]);
% subplot(1,2,2);
% [cs,h]=contourf(X',Y',amp_max,V);
% hold on
% set(gca, 'CLim', [0,1]);
% colorbar('SouthOutside');
% set(gca,'FontSize',18);
% xlabel('k_s \Delta t','FontSize',18); 
% ylabel('k_f \Delta t','FontSize',18); 
% set(gca,'XTick',[ksdtmin:1:ksdtmax]);
% set(gca,'YTick',[kfdtmin:1:kfdtmax]);
% axis([ksdtmin ksdtmax kfdtmin kfdtmax]);
% title(text);