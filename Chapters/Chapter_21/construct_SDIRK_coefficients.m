function [alpha,beta,stages] = construct_SDIRK_coefficients(ti_method)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%IRK2
if (ti_method == 1)
    stages=2;
    alpha=zeros(stages,stages);
    beta=zeros(stages,1);
    alpha(2,1)=0;
    alpha(2,2)=1;
    beta(:)=alpha(stages,:);
elseif (ti_method == 2)
    stages=3;
    alpha=zeros(stages,stages);
    beta=zeros(stages,1);
    alpha(2,1)=1 - 1/sqrt(2);
    alpha(2,2)=alpha(2,1);
    alpha(3,1)=1/(2*sqrt(2));
    alpha(3,2)=alpha(3,1);
    alpha(3,3)=alpha(2,2);
    beta(:)=alpha(stages,:);
elseif (ti_method == 3)
    stages=4;
    alpha=zeros(stages,stages);
    beta=zeros(stages,1);
    alpha(2,1)=1767732205903.0/4055673282236.0;
    alpha(2,2)=1767732205903.0/4055673282236.0;
    alpha(3,1)=2746238789719.0/10658868560708.0;
    alpha(3,2)=-640167445237.0/6845629431997.0;
    alpha(3,3)=alpha(2,2);
    alpha(4,1)=1471266399579.0/7840856788654.0;
    alpha(4,2)=-4482444167858.0/7529755066697.0;
    alpha(4,3)=11266239266428.0/11593286722821.0;
    alpha(4,4)=alpha(2,2);
    beta(:)=alpha(stages,:);
elseif (ti_method == 4)
    stages=6;
    alpha=zeros(stages,stages);
    beta=zeros(stages,1);
    alpha(2,1)=1.0/4.0;
    alpha(2,2)=1.0/4.0;
    alpha(3,1)=8611.0/62500.0;
    alpha(3,2)=-1743.0/31250.0;
    alpha(3,3)=alpha(2,2);
    alpha(4,1)=5012029.0/34652500.0;
    alpha(4,2)=-654441.0/2922500.0;
    alpha(4,3)=174375.0/388108.0;
    alpha(4,4)=alpha(2,2);
    alpha(5,1)=15267082809.0/155376265600.0;
    alpha(5,2)=-71443401.0/120774400.0;
    alpha(5,3)=730878875.0/902184768.0;
    alpha(5,4)=2285395.0/8070912.0;
    alpha(5,5)=alpha(2,2);
    alpha(6,1)=82889.0/524892.0;
    alpha(6,2)=0.0;
    alpha(6,3)=15625.0/83664.0;
    alpha(6,4)=69875.0/102672.0;
    alpha(6,5)=-2260.0/8211.0;
    alpha(6,6)=alpha(2,2);
    beta(:)=alpha(stages,:);
elseif (ti_method == 5)
    stages=8;
    alpha=zeros(stages,stages);
    beta=zeros(stages,1);
    alpha(2,1)=41.0/200.0;
    alpha(2,2)=41.0/200.0;
    alpha(3,1)=41.0/400.0;
    alpha(3,2)=-567603406766.0/11931857230679.0;
    alpha(3,3)=alpha(2,2);
    alpha(4,1)=683785636431.0/9252920307686.0;
    alpha(4,2)=0.0;
    alpha(4,3)=-110385047103.0/1367015193373.0;
    alpha(4,4)=alpha(2,2);
    alpha(5,1)= 3016520224154.0/10081342136671.0;
    alpha(5,2)=0.0;
    alpha(5,3)=30586259806659.0/12414158314087.0;
    alpha(5,4)=-22760509404356.0/11113319521817.0;
    alpha(5,5)=alpha(2,2);
    alpha(6,1)=218866479029.0/1489978393911.0;
    alpha(6,2)=0.0;
    alpha(6,3)=638256894668.0/5436446318841.0;
    alpha(6,4)=-1179710474555.0/5321154724896.0;
    alpha(6,5)=-60928119172.0/8023461067671.0;
    alpha(6,6)=alpha(2,2);
    alpha(7,1)=1020004230633.0/5715676835656.0;
    alpha(7,2)=0.0;
    alpha(7,3)=25762820946817.0/25263940353407.0;
    alpha(7,4)=-2161375909145.0/9755907335909.0;
    alpha(7,5)=-211217309593.0/5846859502534.0;
    alpha(7,6)=-4269925059573.0/7827059040749.0;
    alpha(7,7)=alpha(2,2);
    alpha(8,1)=-872700587467.0/9133579230613.0;
    alpha(8,2)=0.0;
    alpha(8,3)=0.0;
    alpha(8,4)=22348218063261.0/9555858737531.0;
    alpha(8,5)=-1143369518992.0/8141816002931.0;
    alpha(8,6)=-39379526789629.0/19018526304540.0;
    alpha(8,7)=32727382324388.0/42900044865799.0;
    alpha(8,8)=alpha(2,2);
    beta(:)=alpha(stages,:);
elseif (ti_method == 6)
    stages=10;
    alpha=zeros(stages,stages);
    beta=zeros(stages,1);
    alpha(2,2)=0.2928932188134525;
    alpha(3,2)=0.3812815664617709;
    alpha(4,2)=0.4696699141100893;
    alpha(5,2)=0.5580582617584078;
    alpha(6,2)=0.6464466094067263;
    alpha(7,2)=0.7348349570550446;
    alpha(8,2)=0.8232233047033631;
    alpha(9,2)=0.9116116523516815;
    alpha(10,2)=0.7071067811865476;
    alpha(10,10)=0.2928932188134525;
    beta(:)=alpha(stages,:);
end

end

