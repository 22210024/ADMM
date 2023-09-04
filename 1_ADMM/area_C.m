function [Xc,Ploss,Volta,Qg] = area_C(lagrant_c,xigma,Xc_last)
%AREA_B ÇøÓòb,½Úµã5,25-32
%% V=12.66kV S=1000kVA
area = [5,25,26,27,28,29,30,31,32];
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
Xc = [P(25);Q(25);U(5);U(25)];
C = [];
C = [C,P(26) == P(25)-Pload(25)+Pg(25),
    P(27) == P(26)-Pload(26)+Pg(26),
    P(28) == P(27)-Pload(27)+Pg(27),
    P(29) == P(28)-Pload(28)+Pg(28),
    P(30) == P(29)-Pload(29)+Pg(29),
    P(31) == P(30)-Pload(30)+Pg(30),
    P(32) == P(31)-Pload(31)+Pg(31),
    0 == P(32)-Pload(32)+Pg(32)];
C = [C,Q(26) == Q(25)-Qload(25)+Qg(25),
    Q(27) == Q(26)-Qload(26)+Qg(26),
    Q(28) == Q(27)-Qload(27)+Qg(27),
    Q(29) == Q(28)-Qload(28)+Qg(28),
    Q(30) == Q(29)-Qload(29)+Qg(29),
    Q(31) == Q(30)-Qload(30)+Qg(30),
    Q(32) == Q(31)-Qload(31)+Qg(31),
    0 == Q(32)-Qload(32)+Qg(32)];
C = [C,U(33) == 1,
    U(25) == U(5)-2*(r(25)*P(25)+x(25)*Q(25)),
    U(26) == U(25)-2*(r(26)*P(26)+x(26)*Q(26)),
    U(27) == U(26)-2*(r(27)*P(27)+x(27)*Q(27)),
    U(28) == U(27)-2*(r(28)*P(28)+x(28)*Q(28)),
    U(29) == U(28)-2*(r(29)*P(29)+x(29)*Q(29)),
    U(30) == U(29)-2*(r(30)*P(30)+x(30)*Q(30)),
    U(31) == U(30)-2*(r(31)*P(31)+x(31)*Q(31)),
    U(32) == U(31)-2*(r(32)*P(32)+x(32)*Q(32))];
C = [C,-Qgm <= Qg,
    Qg <= Qgm,
    0.95^2 <= U <= 1.05^2];
%% solve
Ploss = r.*(P.^2 + Q.^2);
Objective = sum(Ploss(area(2:end))) + xigma/2*norm(Xc-Xc_last+lagrant_c,2)^2;
obs = sdpsettings('solver','cplex');
optimize(C,Objective,obs)
Ploss = double(Ploss) * S;
P = S*double(P);
Q = S*double(Q);
Qg = S*double(Qg);
U = double(U);
Volta = sqrt(U);
Xc = double(Xc);
end

