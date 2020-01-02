%---------------------------------------------------------------------%
%This code plots the Stability of Time-Differencing Schemes.
%Written by F.X. Giraldo on 10/07
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey, CA 93943-5216
%Variables:
%---------------------------------------------------------------------%

clear; clear all; clf; 
close all;

%Read Number of Files
iplot=0;
% delta=input(' 0=Explicit \n 1=Implicit \n Enter delta: ');
%method=input(' Enter SDIRK Method: ');
method=input(' 1=ARK(1,2,1) \n 2=ARKA(2,3,2) \n 21=ARKB(2,3,2) \n 3=ARK(3,4,3) \n 4=ARK(4,6,4) \n 5=ARK(5,8,5)\n 6=ARK(6,10,6) \n Enter Method: ');
delta=input(' Enter delta (=0 explicit and =1 IMEX): ');
% delta=1;
   
%setup TDS coefficients
text=['ARK' num2str(method)];

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
amp=zeros(ns,nf);
kfdt=zeros(nf,1);
ksdt=zeros(ns,1);

if (method == 0)
    ARK_s=2;
elseif (method == 1)
    ARK_s=2;
elseif (method == 2)
    ARK_s=3;     
elseif (method == 21)
    ARK_s=3;     
elseif (method == 3)
    ARK_s=4;   
elseif (method == 4)
    ARK_s=6;  
elseif (method == 5)
    ARK_s=8; 
elseif (method == 6)
    ARK_s=10; %BV24 Method
end

ARK_A=zeros(ARK_s,ARK_s);
ARK_At=zeros(ARK_s,ARK_s);
ARK_b=zeros(ARK_s,1);
ARK_bt=zeros(ARK_s,1);

if (method == 0)
    %ARK(1,2,1)
    ARK_A(2,1)=1;
    ARK_At(2,2)=1;
    ARK_b(1:ARK_s)= [0,1];
    ARK_bt(1:ARK_s)= [0,1];
elseif (method == 1)
    %ARK(1,2,1) is ARK(1,1,1) but modified to satisfy b=bt.
    ARK_A(1,1)=0;
    ARK_A(1,2)=0;
    ARK_A(2,1)=1;
    ARK_A(2,2)=0;
    ARK_At(2,1)=0;
    ARK_At(2,2)=1;
    ARK_b(1:ARK_s)= ARK_At(ARK_s,1:ARK_s);
    ARK_bt(1:ARK_s)= ARK_b(1:ARK_s);
elseif (method == 2) 
    %ARK2A(2,3,2)
    ARK_A(1,1)=0;
    ARK_A(1,2)=0;
    ARK_A(1,3)=0;
    ARK_A(2,1)=0.5857864376269049511983112757903019214303281246230519;
    ARK_A(2,2)=0;
    ARK_A(2,3)=0;
    ARK_A(3,1)=0.02859547920896831706610375859676730714344270820768398;
    ARK_A(3,2)=0.9714045207910316829338962414032326928565572917923160;
    ARK_A(3,3)=0;
    ARK_At(2,1)=0.292893218813452475599155637895150960715164062311526;
    ARK_At(2,2)=0.2928932188134524755991556378951509607151640623115260;
    ARK_At(3,1)=0.3535533905932737622004221810524245196424179688442370;
    ARK_At(3,2)=0.3535533905932737622004221810524245196424179688442370;
    ARK_At(3,3)=0.2928932188134524755991556378951509607151640623115260;
    ARK_b(1:ARK_s)= ARK_At(ARK_s,1:ARK_s);
    ARK_bt(1:ARK_s)= ARK_b(1:ARK_s);
elseif (method == 21)
    %ARK2B(2,3,2)
    ARK_A(1,1)=0;
    ARK_A(1,2)=0;
    ARK_A(1,3)=0;
    ARK_A(2,1)=0.5857864376269049511983112757903019214303281246230519;
    ARK_A(2,2)=0;
    ARK_A(2,3)=0;
    ARK_A(3,1)=0.5;
    ARK_A(3,2)=0.5;
    ARK_A(3,3)=0;
    ARK_At(2,1)=0.292893218813452475599155637895150960715164062311526;
    ARK_At(2,2)=0.2928932188134524755991556378951509607151640623115260;
    ARK_At(3,1)=0.3535533905932737622004221810524245196424179688442370;
    ARK_At(3,2)=0.3535533905932737622004221810524245196424179688442370;
    ARK_At(3,3)=0.2928932188134524755991556378951509607151640623115260;
    ARK_b(1:ARK_s)= ARK_At(ARK_s,1:ARK_s);
    ARK_bt(1:ARK_s)= ARK_b(1:ARK_s);
elseif (method == 3)  
    %ARK(3,4,3)
    ARK_A(2,1)=1767732205903.0/2027836641118.0;
    ARK_A(3,1)=5535828885825.0/10492691773637.0;
    ARK_A(3,2)=788022342437.0/10882634858940.0;
    ARK_A(4,1)=6485989280629.0/16251701735622.0;
    ARK_A(4,2)=-4246266847089.0/9704473918619.0;
    ARK_A(4,3)=10755448449292.0/10357097424841.0;   
    ARK_At(2,1)=1767732205903.0/4055673282236.0;
    ARK_At(2,2)=1767732205903.0/4055673282236.0;
    ARK_At(3,1)=2746238789719.0/10658868560708.0;
    ARK_At(3,2)=-640167445237.0/6845629431997.0;
    ARK_At(3,3)=ARK_At(2,2);
    ARK_At(4,1)=1471266399579.0/7840856788654.0;
    ARK_At(4,2)=-4482444167858.0/7529755066697.0;
    ARK_At(4,3)=11266239266428.0/11593286722821.0;
    ARK_At(4,4)=ARK_At(2,2);
    ARK_b(1:ARK_s)= ARK_At(ARK_s,1:ARK_s);
    ARK_bt(1:ARK_s)= ARK_b(1:ARK_s);
elseif (method == 4)   
    %ARK(4,6,4)
    ARK_A(2,1)=1.0/2.0;
    ARK_A(3,1)=13861.0/62500.0;
    ARK_A(3,2)=6889.0/62500.0;
    ARK_A(4,1)=-116923316275.0/2393684061468.0;
    ARK_A(4,2)=-2731218467317.0/15368042101831.0;
    ARK_A(4,3)=9408046702089.0/11113171139209.0;
    ARK_A(5,1)=-451086348788.0/2902428689909.0;
    ARK_A(5,2)=-2682348792572.0/7519795681897.0;
    ARK_A(5,3)=12662868775082.0/11960479115383.0;
    ARK_A(5,4)=3355817975965.0/11060851509271.0;
    ARK_A(6,1)=647845179188.0/3216320057751.0;
    ARK_A(6,2)=73281519250.0/8382639484533.0;
    ARK_A(6,3)=552539513391.0/3454668386233.0;
    ARK_A(6,4)=3354512671639.0/8306763924573.0;
    ARK_A(6,5)=4040.0/17871.0;
    ARK_At(2,1)=1.0/4.0;
    ARK_At(2,2)=1.0/4.0;
    ARK_At(3,1)=8611.0/62500.0;
    ARK_At(3,2)=-1743.0/31250.0;
    ARK_At(3,3)=ARK_At(2,2);
    ARK_At(4,1)=5012029.0/34652500.0;
    ARK_At(4,2)=-654441.0/2922500.0;
    ARK_At(4,3)=174375.0/388108.0;
    ARK_At(4,4)=ARK_At(2,2);
    ARK_At(5,1)=15267082809.0/155376265600.0;
    ARK_At(5,2)=-71443401.0/120774400.0;
    ARK_At(5,3)=730878875.0/902184768.0;
    ARK_At(5,4)=2285395.0/8070912.0;
    ARK_At(5,5)=ARK_At(2,2);
    ARK_At(6,1)=82889.0/524892.0;
    ARK_At(6,2)=0.0;
    ARK_At(6,3)=15625.0/83664.0;
    ARK_At(6,4)=69875.0/102672.0;
    ARK_At(6,5)=-2260.0/8211.0;
    ARK_At(6,6)=ARK_At(2,2);
    ARK_b(1:ARK_s)= ARK_At(ARK_s,1:ARK_s);
    ARK_bt(1:ARK_s)=ARK_b(1:ARK_s);    
elseif (method == 5) 
    %ARK(5,8,5)
    ARK_A(2,1)=41.0/100.0;
    ARK_A(3,1)=367902744464.0/2072280473677.0;
    ARK_A(3,2)=677623207551.0/8224143866563.0;
    ARK_A(4,1)=1268023523408.0/10340822734521.0;
    ARK_A(4,2)=0.0;
    ARK_A(4,3)=1029933939417.0/13636558850479.0;
    ARK_A(5,1)=14463281900351.0/6315353703477.0;
    ARK_A(5,2)=0.0;
    ARK_A(5,3)=66114435211212.0/5879490589093.0;
    ARK_A(5,4)=-54053170152839.0/4284798021562.0;
    ARK_A(6,1)=14090043504691.0/34967701212078.0;
    ARK_A(6,2)=0.0;
    ARK_A(6,3)=15191511035443.0/11219624916014.0;
    ARK_A(6,4)=-18461159152457.0/12425892160975.0;
    ARK_A(6,5)=-281667163811.0/9011619295870.0;
    ARK_A(7,1)=19230459214898.0/13134317526959.0;
    ARK_A(7,2)=0.0;
    ARK_A(7,3)=21275331358303.0/2942455364971.0;
    ARK_A(7,4)=-38145345988419.0/4862620318723.0;
    ARK_A(7,5)=-1.0/8.0;
    ARK_A(7,6)=-1.0/8.0;
    ARK_A(8,1)=-19977161125411.0/11928030595625.0;
    ARK_A(8,2)=0.0;
    ARK_A(8,3)=-40795976796054.0/6384907823539.0;
    ARK_A(8,4)=177454434618887.0/12078138498510.0;
    ARK_A(8,5)=782672205425.0/8267701900261.0;
    ARK_A(8,6)=-69563011059811.0/9646580694205.0;
    ARK_A(8,7)=7356628210526.0/4942186776405.0;
    ARK_At(2,1)=41.0/200.0;
    ARK_At(2,2)=41.0/200.0;
    ARK_At(3,1)=41.0/400.0;
    ARK_At(3,2)=-567603406766.0/11931857230679.0;
    ARK_At(3,3)=ARK_At(2,2);
    ARK_At(4,1)=683785636431.0/9252920307686.0;
    ARK_At(4,2)=0.0;
    ARK_At(4,3)=-110385047103.0/1367015193373.0;
    ARK_At(4,4)=ARK_At(2,2);
    ARK_At(5,1)= 3016520224154.0/10081342136671.0;
    ARK_At(5,2)=0.0;
    ARK_At(5,3)=30586259806659.0/12414158314087.0;
    ARK_At(5,4)=-22760509404356.0/11113319521817.0;
    ARK_At(5,5)=ARK_At(2,2);
    ARK_At(6,1)=218866479029.0/1489978393911.0;
    ARK_At(6,2)=0.0;
    ARK_At(6,3)=638256894668.0/5436446318841.0;
    ARK_At(6,4)=-1179710474555.0/5321154724896.0;
    ARK_At(6,5)=-60928119172.0/8023461067671.0;
    ARK_At(6,6)=ARK_At(2,2);
    ARK_At(7,1)=1020004230633.0/5715676835656.0;
    ARK_At(7,2)=0.0;
    ARK_At(7,3)=25762820946817.0/25263940353407.0;
    ARK_At(7,4)=-2161375909145.0/9755907335909.0;
    ARK_At(7,5)=-211217309593.0/5846859502534.0;
    ARK_At(7,6)=-4269925059573.0/7827059040749.0;
    ARK_At(7,7)=ARK_At(2,2);
    ARK_At(8,1)=-872700587467.0/9133579230613.0;
    ARK_At(8,2)=0.0;
    ARK_At(8,3)=0.0;
    ARK_At(8,4)=22348218063261.0/9555858737531.0;
    ARK_At(8,5)=-1143369518992.0/8141816002931.0;
    ARK_At(8,6)=-39379526789629.0/19018526304540.0;
    ARK_At(8,7)=32727382324388.0/42900044865799.0;
    ARK_At(8,8)=ARK_At(2,2);
    ARK_b(1:ARK_s)= ARK_At(ARK_s,1:ARK_s);
    ARK_bt(1:ARK_s)=ARK_b(1:ARK_s);
elseif (method == 6)  
    %ARK(6,10,6)
    ARK_A(2,1)=0.2928932188134525;
    ARK_A(3,1)=0.3602847895715037;
    ARK_A(3,2)=0.02099677689026724;
    ARK_A(4,1)=0.4267095308101442;
    ARK_A(4,3)=0.04296038329994516;
    ARK_A(5,1)=0.4901881062634202;
    ARK_A(5,4)=0.06787015549498755;
    ARK_A(6,1)=0.548902422602367;
    ARK_A(6,5)=0.09754418680435926;
    ARK_A(7,1)=0.6003576019251088;
    ARK_A(7,6)=0.1344773551299359;
    ARK_A(8,1)=0.6404229388942057;
    ARK_A(8,7)=0.1828003658091574;
    ARK_A(9,1)=0.6612274743186617;
    ARK_A(9,8)=0.2503841780330199;    
    ARK_A(10,1)=0.3942428457889902;
    ARK_A(10,9)=0.6057571542110098;
    ARK_At(2,2)=0.2928932188134525;
    ARK_At(3,2)=0.3812815664617709;
    ARK_At(4,2)=0.4696699141100893;
    ARK_At(5,2)=0.5580582617584078;
    ARK_At(6,2)=0.6464466094067263;
    ARK_At(7,2)=0.7348349570550446;
    ARK_At(8,2)=0.8232233047033631;
    ARK_At(9,2)=0.9116116523516815;
    ARK_At(10,2)=0.7071067811865476;
    ARK_At(10,10)=0.2928932188134525;
    ARK_b(1:ARK_s)= ARK_At(ARK_s,1:ARK_s);
    ARK_bt(1:ARK_s)=ARK_b(1:ARK_s);
end

for l=1:nf
   kfdt(l)=kfdtmin + (l-1)*dkfdt;
   for k=1:ns
      ksdt(k)=ksdtmin + (k-1)*dksdt;
      q=ones(ARK_s,1);
      for i=2:ARK_s
          %Explicit Part
          for j=1:i-1
            q(i)=q(i) + 1i*( ARK_A(i,j)*ksdt(k) + ARK_At(i,j)*kfdt(l) )*q(j);
          end
          %Implicit Part
          q(i)=q(i)/(1 - delta*1i*kfdt(l)*ARK_At(i,i));
      end
      temp=1;
      for i=1:ARK_s
          temp=temp + 1i*ARK_b(i)*q(i)*(ksdt(k) + kfdt(l));
      end
      amp1d(k,l)=sqrt( temp*conj(temp) );
      amp(k,l)=amp1d(k,l);
      if (amp(k,l) > 1 + 1e-6) 
          amp(k,l)=nan;
      end
   end       
end

[X,Y]=meshgrid(ksdt,kfdt);
V=0:0.05:1.05;
amp1=zeros(ns,nf);
amp2=zeros(ns,nf);
amp_max=zeros(ns,nf);
% figure
% subplot(1,2,1);
% plot(ksdt,amp1d(:,1),'r-','LineWidth',2);
% set(gca,'FontSize',18);
% xlabel('k_s \Delta t','FontSize',18); 
% ylabel('|A|','FontSize',18); 
% %axis([-cmax cmax 0 1.5]);
% axis([-1 1 0 1.5]);
% subplot(1,2,2);
% [cs,h]=contourf(X',Y',amp,V);
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

% figure
% subplot(1,2,1);
% plot(ksdt,amp(:,1),'r-','LineWidth',2);
% set(gca,'FontSize',18);
% xlabel('k_s \Delta t','FontSize',18); 
% ylabel('|A|','FontSize',18); 
% axis([ksdtmin ksdtmax 0 1.75]);
% subplot(1,2,2);
% colormap(jet);
% [cs,h]=contourf(X',Y',amp,V);
% hold on
% set(gca, 'CLim', [0,1]);
% colorbar('SouthOutside');
% set(gca,'FontSize',18);
% xlabel('k_s \Delta t','FontSize',18); 
% ylabel('k_f \Delta t','FontSize',18); 
% set(gca,'XTick',[ksdtmin:1:ksdtmax]);
% set(gca,'YTick',[kfdtmin:10:kfdtmax]);
% axis([ksdtmin ksdtmax kfdtmin kfdtmax]);
% title(text);

figure
colormap(jet);
[cs,h]=contourf(X',Y',amp,V);
hold on
set(gca, 'CLim', [0,1]);
colorbar('SouthOutside');
set(gca,'FontSize',18);
xlabel('k_s \Delta t','FontSize',18); 
ylabel('k_f \Delta t','FontSize',18); 
set(gca,'XTick',[ksdtmin:1:ksdtmax]);
set(gca,'YTick',[kfdtmin:10:kfdtmax]);
axis([ksdtmin ksdtmax kfdtmin kfdtmax]);
title(text);