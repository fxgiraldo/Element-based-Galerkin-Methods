%---------------------------------------------------------------------%
%This function computes the LSRK Method by Carpenter-Kennedy 1994.
%Written by F.X. Giraldo on April 19, 2019
%           Department of Applied Mathematics
%           Naval Postgraduate School
%           Monterey; CA 93943-5216
%---------------------------------------------------------------------%
function [q0] = ti_LSRK(q0,u,Dhat,Fhat,intma,periodicity,time,ntime,dt)

%Initialize RK coefficients
RKA = [(0), 
       (-567301805773) / (1357537059087), 
       (-2404267990393) / (2016746695238), 
       (-3550918686646) / (2091501179385), 
       (-1275806237668) / (842570457699 )];

RKB = [(1432997174477) / (9575080441755 ),
       (5161836677717) / (13612068292357),
       (1720146321549) / (2090206949498 ),
       (3134564353537) / (4481467310338 ),
       (2277821191437) / (14882151754819)];

RKC = [(0),
       (1432997174477) / (9575080441755),
       (2526269341429) / (6820363962896),
       (2006345519317) / (3224310063776),
       (2802321613138) / (2924317926251)];

Npoin=length(q0);
dq=zeros(Npoin,1);
qp=q0;
stages=length(RKA);

%Time Integration
for itime=1:ntime
    time=time + dt;
    %RK Stages
    for s = 1:stages
        %Create RHS Matrix
        %R = create_rhs(qp,u,Dhat,Fhat,intma);
        R=Dhat*qp*u;
        %Solve System
        for I=1:Npoin
            dq(I) = RKA(s)*dq(I) + dt*R(I);
            qp(I) = qp(I) + RKB(s)*dq(I);
        end
        if periodicity(Npoin) == periodicity(1)
            qp(Npoin)=qp(1); %periodicity
        end
    end %s
    %Update Q
    q0 = qp;
end %itime
