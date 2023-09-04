function [Xb,Ploss,Volta,Qg] = area_B(lagrant_b,xigma,Xb_last)
%AREA_B ÇøÓòb,½Úµã6-17
%% V=12.66kV S=1000kVA
area = [6,7,8,9,10,11,12,13,14,15,16,17];
mpc = case33bw;
r = mpc.branch(:,3);
x = mpc.branch(:,4);
Pload = mpc.bus(2:end,3);
Qload = mpc.bus(2:end,4);
Pg = zeros(33,1);
genbus = [7,11,15,20,29];
Pg(genbus) = 400;
Qgm = zeros(33,1);
Qgm(genbus) = 400; 
V = 12.66;S = 1000;r = r/V^2; x = x/V^2; Pload = Pload/S; Qload = Qload/S; Pg = Pg/S; Qgm = Qgm/S;
%% variables and constraints
P = sdpvar(32,1);
Q = sdpvar(32,1);
Qg = sdpvar(33,1);
U = sdpvar(33,1);
Xb = [P(7);Q(7);U(6);U(7)];
C = [];
C = [C,P(8) == P(7)-Pload(7)+Pg(7),
    P(9) == P(8)-Pload(8)+Pg(8),
    P(10) == P(9)-Pload(9)+Pg(9),
    P(11) == P(10)-Pload(10)+Pg(10),
    P(12) == P(11)-Pload(11)+Pg(11),
    P(13) == P(12)-Pload(12)+Pg(12),
    P(14) == P(13)-Pload(13)+Pg(13),
    P(15) == P(14)-Pload(14)+Pg(14),
    P(16) == P(15)-Pload(15)+Pg(15),
    P(17) == P(16)-Pload(16)+Pg(16),
    0 == P(17)-Pload(17)+Pg(17)];
C = [C,Q(8) == Q(7)-Qload(7)+Qg(7),
    Q(9) == Q(8)-Qload(8)+Qg(8),
    Q(10) == Q(9)-Qload(9)+Qg(9),
    Q(11) == Q(10)-Qload(10)+Qg(10),
    Q(12) == Q(11)-Qload(11)+Qg(11),
    Q(13) == Q(12)-Qload(12)+Qg(12),
    Q(14) == Q(13)-Qload(13)+Qg(13),
    Q(15) == Q(14)-Qload(14)+Qg(14),
    Q(16) == Q(15)-Qload(15)+Qg(15),
    Q(17) == Q(16)-Qload(16)+Qg(16),
    0 == Q(17)-Qload(17)+Qg(17)];
C = [C,U(33) == 1,
    U(7) == U(6)-2*(r(7)*P(7)+x(7)*Q(7)),
    U(8) == U(7)-2*(r(8)*P(8)+x(8)*Q(8)),
    U(9) == U(8)-2*(r(9)*P(9)+x(9)*Q(9)),
    U(10) == U(9)-2*(r(10)*P(10)+x(10)*Q(10)),
    U(11) == U(10)-2*(r(11)*P(11)+x(11)*Q(11)),
    U(12) == U(11)-2*(r(12)*P(12)+x(12)*Q(12)),
    U(13) == U(12)-2*(r(13)*P(13)+x(13)*Q(13)),
    U(14) == U(13)-2*(r(14)*P(14)+x(14)*Q(14)),
    U(15) == U(14)-2*(r(15)*P(15)+x(15)*Q(15)),
    U(16) == U(15)-2*(r(16)*P(16)+x(16)*Q(16)),
    U(17) == U(16)-2*(r(17)*P(17)+x(17)*Q(17))];
C = [C,-Qgm <= Qg,
    Qg <= Qgm,
    0.95^2 <= U <= 1.05^2];
%% solve
Ploss = r.*(P.^2 + Q.^2);
Objective = sum(Ploss(area(2:end))) + xigma/2*norm(Xb-Xb_last+lagrant_b,2)^2;
obs = sdpsettings('solver','cplex');
optimize(C,Objective,obs);
Ploss = double(Ploss) * S;
P = S*double(P);
Q = S*double(Q);
Qg = S*double(Qg);
U = double(U);
Volta = sqrt(U);
Xb = double(Xb);
end

