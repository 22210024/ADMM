function [Xa,Ploss,Volta,Qg]=area_A(lagrant_a,xigma,Xa_last)
%AREA_A 区域节点33,1-7,18-25
%% V=12.66kV S=1000kVA
area = [33,1,2,3,4,5,6,7,18,19,20,21,22,23,24,25];%计算Ploss时首尾不考虑
mpc = case33bw;
r = mpc.branch(:,3);
x = mpc.branch(:,4);
Pload = mpc.bus(2:end,3);
Qload= mpc.bus(2:end,4);
Pg = zeros(32,1);
genbus=[7,11,15,20,29];
Pg(genbus) = 400;
Qgm = zeros(33,1);
Qgm(genbus)= 400;
V = 12.66;S = 1000;r = r/V^2;x = x/V^2;Pload = Pload/S; Qload = Qload/S;Pg = Pg/S;Qgm = Qgm/S;
%% variable and constraints
P = sdpvar(32,1);
Q = sdpvar(32,1);
Qg = sdpvar(33,1);
U = sdpvar(33,1);
Xa = [P(7);Q(7);U(6);U(7);P(25);Q(25);U(5);U(25)];
C = [];
C = [C,P(2)+P(18) == P(1)-Pload(1)+Pg(1),
    P(3)+P(22) == P(2)-Pload(2)+Pg(2),
    P(4) == P(3)-Pload(3)+Pg(3),
    P(5) == P(4)-Pload(4)+Pg(4),
    P(6)+P(25) == P(5)-Pload(5)+Pg(5),
    P(7) == P(6)-Pload(6)+Pg(6),
    P(19) == P(18)-Pload(18)+Pg(18),
    P(20) == P(19)-Pload(19)+Pg(19),
    P(21) == P(20)-Pload(20)+Pg(20),
    0 == P(21)-Pload(21)+Pg(21),
    P(23) == P(22)-Pload(22)+Pg(22),
    P(24) == P(23)-Pload(23)+Pg(23),
    0 == P(24)-Pload(24)+Pg(24)];
C = [C,Q(2)+Q(18) == Q(1)-Qload(1)+Qg(1),
    Q(3)+Q(22) == Q(2)-Qload(2)+Qg(2),
    Q(4) == Q(3)-Qload(3)+Qg(3),
    Q(5) == Q(4)-Qload(4)+Qg(4),
    Q(6)+Q(25) == Q(5)-Qload(5)+Qg(5),
    Q(7) == Q(6)-Qload(6)+Qg(6),
    Q(19) == Q(18)-Qload(18)+Qg(18),
    Q(20) == Q(19)-Qload(19)+Qg(19),
    Q(21) == Q(20)-Qload(20)+Qg(20),
    0 == Q(21)-Qload(21)+Qg(21),
    Q(23) == Q(22)-Qload(22)+Qg(22),
    Q(24) == Q(23)-Qload(23)+Qg(23),
    0 == Q(24)-Qload(24)+Qg(24)];
C = [C,U(33) == 1,
    U(1) == U(33)-2*(r(1)*P(1)+x(1)*Q(1)),
    U(2) == U(1)-2*(r(2)*P(2)+x(2)*Q(2)),
    U(3) == U(2)-2*(r(3)*P(3)+x(3)*Q(3)),
    U(4) == U(3)-2*(r(4)*P(4)+x(4)*Q(4)),
    U(5) == U(4)-2*(r(5)*P(5)+x(5)*Q(5)),
    U(6) == U(5)-2*(r(6)*P(6)+x(6)*Q(6)),
    U(7) == U(6)-2*(r(7)*P(7)+x(7)*Q(7)),
    U(18) == U(1)-2*(r(18)*P(18)+x(18)*Q(18)),
    U(19) == U(18)-2*(r(19)*P(19)+x(19)*Q(19)),
    U(20) == U(19)-2*(r(20)*P(20)+x(20)*Q(20)),
    U(21) == U(20)-2*(r(21)*P(21)+x(21)*Q(21)),
    U(22) == U(2)-2*(r(22)*P(22)+x(22)*Q(22)),
    U(23) == U(22)-2*(r(23)*P(23)+x(23)*Q(23)),
    U(24) == U(23)-2*(r(24)*P(24)+x(24)*Q(24))];
C = [C,-Qgm <= Qg,
    Qg <= Qgm,
    0.95^2<= U <=1.05^2];
%% solve
Ploss = r.*(P.^2+Q.^2);
Objective = sum(Ploss(area(2:15))) + xigma/2*norm(Xa - Xa_last + lagrant_a, 2)^2;%norm(Xa - Xa_last + lagrant_a, 2);two norm 
ops = sdpsettings('solver', 'cplex');
optimize(C, Objective, ops) % Set the solver to cplex
Ploss = double(Ploss) * S;
P = S*double(P);
Q = S*double(Q);
Qg = S*double(Qg);
U = double(U);
Volta = sqrt(U);
Xa=double(Xa);
end

