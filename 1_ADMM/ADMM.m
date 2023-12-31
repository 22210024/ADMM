lagrant_a = [0;0;0;0;0;0;0;0];xigma = 0;Xa_last = [0;0;0;0;0;0;0;0];
[Xa,Ploss_a,Volta_a,Qg_a]=area_A(lagrant_a,xigma,Xa_last);
lagrant_b = [0;0;0;0];xigma = 0;Xb_last = [0;0;0;0];
[Xb,Ploss_b,Volta_b,Qg_b]=area_B(lagrant_b,xigma,Xb_last);
lagrant_c = [0;0;0;0];xigma = 0;Xc_last = [0;0;0;0];
[Xc,Ploss_c,Volta_c,Qg_c]=area_C(lagrant_c,xigma,Xc_last);
Xbc = [Xb;Xc];
Xa_last = (Xa + Xbc)/2;
Xb_last = Xa_last(1:4);Xc_last = Xa_last(5:8);
lagrant_a = lagrant_a + (Xa -Xa_last);
lagrant_b = lagrant_b + (Xb -Xb_last);
lagrant_c = lagrant_c + (Xc -Xc_last);
xigma = 0.05;
Anode = [33,1,2,3,4,5,6,7,18,19,20,21,22,23,24,25];
Bnode = [6,7,8,9,10,11,12,13,14,15,16,17];
Cnode = [5,25,26,27,28,29,30,31,32];
Ploss = zeros(32,1);
Ploss(Anode(2:15)) = Ploss_a(Anode(2:15));
Ploss(Bnode(2:end)) = Ploss_b(Bnode(2:end));
Ploss(Cnode(2:end)) = Ploss_c(Cnode(2:end));
Ploss = sum(Ploss);
plot(1,Ploss,'r*')
hold on
plot(1,norm(Xa-[Xb;Xc],2)^2);
hold on
for k = 1:100
    if norm(Xa-[Xb;Xc],2)^2 > 1 / (10^4)
        [Xa,Ploss_a,Volta_a,Qg_a]=area_A(lagrant_a,xigma,Xa_last);
        [Xb,Ploss_b,Volta_b,Qg_b]=area_B(lagrant_b,xigma,Xb_last);
        [Xc,Ploss_c,Volta_c,Qg_c]=area_C(lagrant_c,xigma,Xc_last);
        Xbc = [Xb;Xc];
        Xa_last = (Xa + Xbc)/2;
        Xb_last = Xa_last(1:4);Xc_last = Xa_last(5:8);
        lagrant_a = lagrant_a + (Xa -Xa_last);
        lagrant_b = lagrant_b + (Xb -Xb_last);
        lagrant_c = lagrant_c + (Xc -Xc_last);
        Ploss(Anode(2:15)) = Ploss_a(Anode(2:15));
        Ploss(Bnode(2:end)) = Ploss_b(Bnode(2:end));
        Ploss(Cnode(2:end)) = Ploss_c(Cnode(2:end));
        Ploss = sum(Ploss);
        figure(1)
        plot(k+1,Ploss,'r*')
        hold on
        figure(2)
        plot(k+1,norm(Xa-[Xb;Xc],2)^2,'b*');
        hold on
    end
end
Qg = zeros(33,1);
Qg(Anode(1:15)) = Qg_a(Anode(1:15));
Qg(Bnode(2:end)) = Qg_b(Bnode(2:end));
Qg(Cnode(2:end)) = Qg_c(Cnode(2:end));