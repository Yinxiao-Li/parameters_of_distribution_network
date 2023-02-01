function [P1, Pt1, Ca1, PG, QG, Pac1, Qt1, SOC1, delta1, delta2, delta3]...
    = opt_oneday(SOC0, W, Tf, ct, kload, Cb)
load distribution_network;
mpc.branch(:,5) = 0;
SB = mpc.baseMVA;

l0 = 7.503001200480192e-04;
l2 = 0.025510204081633;

PV_r = 0.02 / SB; 
Tpv = Tf + 0.03 * W;
PV_ava = W / 1000 .* (1 - 0.0045 * (Tpv - 25)) * PV_r; 


nbat = 2424; 
% Cb = 2.5;
Cr = 0.5; 
Cbat = 2.5*3.3*nbat/1000*1800; % £¤2000/kWh
Tref = 298.15;
Eacal = 20592;
Rg = 8.314;
a0 = 0.384;
k0 = 0.142;
Uaref = 0.123;
F = 96485;
ku = 0.139088471038960;
U0 = 3.207526748442694;

kcalref = 5.535e-04;
kltref = 1.203e-03;
khtref = 4.38e-04;

R0 = 0.0976*exp(-0.0077*Tf) + 2.5e-8*Tf.^2+2e-6*Tf+6e-5;
Cb0 = zeros(24,5);
for k = 1:5
    Cb0(:,k) = Cb(1,k) * (-0.354*exp(-0.031*Tf)+1.157);
end

nBESS = 5;
Pt = sdpvar(24,nBESS);
Qt = sdpvar(24,nBESS);
Pac = sdpvar(24,nBESS);
Pcl = sdpvar(24,nBESS);
SOC = sdpvar(24,nBESS);
DSOC = sdpvar(24,nBESS);
Pr = 2.5*3.3; 

P = sdpvar(59,24);
Q = sdpvar(59,24);
I = sdpvar(59,24);
V = sdpvar(60,24);
nPV = 5;
Pg = sdpvar(nPV,24);
Qg = sdpvar(nPV,24);

% ³£Êý
r = mpc.branch(1:59,3);
x = mpc.branch(1:59,4);
Pz = mpc.bus(1:60,3);
Qz = mpc.bus(1:60,4);
Pz = Pz / SB;
Qz = Qz / SB;

C = [SOC >= 0.2, SOC <= 0.9];
C = [C, SOC(24,5) == 0.2 * ones(1,5)];

for j = 1:24
    C = [C, V(1,j) == 1.03^2];
    C = [C, 0.95^2 <= V(2:60,j) <= 1.05^2];
    C = [C, Pg(1,j) >= 0];
    C = [C, Pg(:,j) <= PV_ava(j)];
    for k = 1:nPV
        C = [C, cone([Pg(k,j); Qg(k,j)] , PV_r)];
    end
    
    for i = [(1:6) (8:11) (13:19) (21:35) (38:59)]
        parent = find(mpc.branch(:,1) == mpc.branch(i,2));
        C = [C, P(i,j) == sum(P(parent,j)) + kload(j)*Pz(mpc.branch(i,2),1) + I(i,j) * r(i)];
        C = [C, Q(i,j) == sum(Q(parent,j)) + kload(j)*Qz(mpc.branch(i,2),1) + I(i,j) * x(i)];
    end
    
    i = 7;
    parent = find(mpc.branch(:,1) == mpc.branch(i,2));
    C = [C, P(i,j) == sum(P(parent,j)) + kload(j)*Pz(mpc.branch(i,2),1) + I(i,j) * r(i) - Pg(1,j) + Pac(j,1)*nbat/SB/1000000];
    C = [C, Q(i,j) == sum(Q(parent,j)) + kload(j)*Qz(mpc.branch(i,2),1) + I(i,j) * x(i) - Qg(1,j) + Qt(j,1)*nbat/SB/1000000];
    
    i = 12;
    parent = find(mpc.branch(:,1) == mpc.branch(i,2));
    C = [C, P(i,j) == sum(P(parent,j)) + kload(j)*Pz(mpc.branch(i,2),1) + I(i,j) * r(i) - Pg(2,j) + Pac(j,2)*nbat/SB/1000000];
    C = [C, Q(i,j) == sum(Q(parent,j)) + kload(j)*Qz(mpc.branch(i,2),1) + I(i,j) * x(i) - Qg(2,j) + Qt(j,2)*nbat/SB/1000000];

    i = 20;
    parent = find(mpc.branch(:,1) == mpc.branch(i,2));
    C = [C, P(i,j) == sum(P(parent,j)) + kload(j)*Pz(mpc.branch(i,2),1) + I(i,j) * r(i) - Pg(3,j) + Pac(j,3)*nbat/SB/1000000];
    C = [C, Q(i,j) == sum(Q(parent,j)) + kload(j)*Qz(mpc.branch(i,2),1) + I(i,j) * x(i) - Qg(3,j) + Qt(j,3)*nbat/SB/1000000];

    i = 36;
    parent = find(mpc.branch(:,1) == mpc.branch(i,2));
    C = [C, P(i,j) == sum(P(parent,j)) + kload(j)*Pz(mpc.branch(i,2),1) + I(i,j) * r(i) - Pg(4,j) + Pac(j,4)*nbat/SB/1000000];
    C = [C, Q(i,j) == sum(Q(parent,j)) + kload(j)*Qz(mpc.branch(i,2),1) + I(i,j) * x(i) - Qg(4,j) + Qt(j,4)*nbat/SB/1000000];

    i = 37;
    parent = find(mpc.branch(:,1) == mpc.branch(i,2));
    C = [C, P(i,j) == sum(P(parent,j)) + kload(j)*Pz(mpc.branch(i,2),1) + I(i,j) * r(i) - Pg(5,j) + Pac(j,5)*nbat/SB/1000000];
    C = [C, Q(i,j) == sum(Q(parent,j)) + kload(j)*Qz(mpc.branch(i,2),1) + I(i,j) * x(i) - Qg(5,j) + Qt(j,5)*nbat/SB/1000000];

    for i = 1:59
        C = [C, V(mpc.branch(i,2),j) == V(mpc.branch(i,1),j) - 2 * (r(i) * P(i,j) + x(i) * Q(i,j)) + (r(i)^2 + x(i)^2) * I(i,j)];
        C = [C, cone([2 * P(i,j); 2 * Q(i,j); I(i,j) - V(mpc.branch(i,1),j) ] , I(i,j) + V(mpc.branch(i,1),j))];
    end
    
end

for k = 1:nBESS

    C = C + [-0.5 <= DSOC <= 0.5];
    
% t = 1
    C = C + [SOC0(k,1) + DSOC(1,k) == SOC(1,k)];
    C = C + [Pt(1,k) >= U0*Cb0(1,k)*DSOC(1,k) + ...
    (0.5*ku*Cb0(1,k)+Cb0(1,k)^2*R0(1,1))*(DSOC(1,k)+ku*Cb0(1,k)/2/(0.5*ku*Cb0(1,k)+Cb0(1,k)^2*R0(1,1))*SOC0(k,1))^2 + ...
    (ku*Cb0(1,k))^2/4/(0.5*ku*Cb0(1,k)+Cb0(1,k)^2*R0(1,1))*(0.04-0.4*SOC0(k,1))];
% t = 2-24
    for i = 2:24
        C = C + [SOC(i-1,k) + DSOC(i,k) == SOC(i,k)];
        C = C + [Pt(i,k) >= U0*Cb0(i,k)*DSOC(i,k) + ...
        (0.5*ku*Cb0(i,k)+Cb0(i,k)^2*R0(i,1))*(DSOC(i,k)+ku*Cb0(i,k)/2/(0.5*ku*Cb0(i,k)+Cb0(i,k)^2*R0(i,1))*SOC(i-1,k))^2 + ...
        (ku*Cb0(i,k))^2/4/(0.5*ku*Cb0(i,k)+Cb0(i,k)^2*R0(i,1))*(0.04-0.4*SOC(i-1,k))];
    end

    for j = 1:24
        C = [C, Pac(j,k) == Pt(j,k) + Pcl(j,k)];
        C = [C, Pcl(j,k) >= l0*Pr + l2/Pr*Pac(j,k)^2];
        C = [C, cone([Pac(j,k); Qt(j,k)] , 2.5*3.3*0.5*1.1 )];
    end
end

x1 = (0.2:0.01:0.9)';
xa = 8.5e-3 + x1*(0.78 - 8.5e-3);
Ua = 0.6379 + 0.5416*exp(-305.5309*xa) + 0.044*tanh(-(xa-0.1958)/0.1088) ...
    - 0.1978*tanh((xa-1.0571)/0.0854) - 0.6875*tanh((xa+0.0117)/0.0529) ...
    - 0.0175*tanh((xa-0.5692)/0.0875);
kcals1 = exp(a0*F/Rg*(Uaref-Ua)/Tref) + k0;
X1 = [ones(length(x1),1) x1];
bcal = regress(kcals1,X1);


y1 = (0:0.02:2.5)';
Ta = Tf + 273.15;
Ty1 = Ta(1) + y1.^2*R0(1)*(9.12+3.273);
kcycht1 = (exp(-32699/8.314*(1./Ty1-1/298.15)));
kcyclt1 = (exp(55546/8.314*(1./Ty1-1/298.15)).*exp(2.64*(y1-2.5)/2.5));
kch = (kcycht1*khtref+kcyclt1*kltref).* abs(y1).^0.5;
kch  = kch.*exp(-Eacal/Rg*(1./Ty1-1/Tref));
fun1 = @(bch,y1)bch*y1;
options = optimoptions('lsqcurvefit','Display','off');
bch = lsqcurvefit(fun1,0,y1,kch,-10000,10000,options);

y2 = (-2.5:0.02:0)';
Ty2 = Ta(1) + y2.^2*R0(1)*(9.12+3.273);
kcycht2 = (exp(-32699/8.314*(1./Ty2-1/298.15)));
kdis = exp(-Eacal/Rg*(1./Ty2-1/Tref)).*kcycht2*khtref.* abs(y2).^0.5;
fun2 = @(bdis,y2)bdis*y2;
bdis = lsqcurvefit(fun2,0,y2,kdis,-10000,10000,options);
[x,y] = meshgrid(0.2:0.1:0.9,(-2.5:0.5:2.5));
T = Ta(1) + y.^2*R0(1)*(9.12+3.273);
kcychT = (exp(-32699/8.314*(1./T-1/298.15)));
khT = kcychT*khtref;
kcyclT = (exp(55546/8.314*(1./T-1/298.15)).*exp(2.64*(y-2.5)/2.5));
klT = kcyclT*kltref;
klT(1:(length(klT)+1)/2,:) = 0;
kcycT = khT+klT;
y3 = (-2.5:0.02:2.5)';
funcyc = @(bcyc,y)(bcyc(1)+bcyc(2)*y).*(y<=0)+(bcyc(1)+bcyc(3)*y).*(y>0);
bcyc = lsqcurvefit(funcyc,[0,bdis,bch],y3,[kdis; kch(2:126,1)],[0,-10000,-10000],...
        [10000,10000,10000],options);

ZM = ((bcal(1)+bcal(2)*x(1:6,:)).*exp(-Eacal/Rg*(1./T(1:6,:)-1/Tref))*kcalref).^2 ...
    + kcycT(1:6,:).^2.*(-y(1:6,:));
ZM0 = ZM + 2*kcalref*(bcal(1)*bcyc(2)*y(1:6,:)+bcal(2)*bcyc(1)*x(1:6,:)+bcal(1)*bcyc(1));
ZM1 = ZM0 + 2*kcalref*(-bcyc(2))*bcal(2)*0.2*abs(y(1:6,:));
ZM2 = ZM0 + 2*kcalref*(-bcyc(2))*bcal(2)*(0.9*abs(y(1:6,:))+2.5*x(1:6,:)-0.9*2.5);

ZM(7:11,:) = ((bcal(1)+bcal(2)*x(7:11,:)).*exp(-Eacal/Rg*(1./T(7:11,:)-1/Tref))*kcalref).^2 ...
    + kcycT(7:11,:).^2.*y(7:11,:);
ZM0(7:11,:) = ZM(7:11,:) + 2*kcalref*(bcal(1)*bcyc(3)*y(7:11,:)+bcal(2)*bcyc(1)*x(7:11,:)+bcal(1)*bcyc(1));
ZM1(7:11,:) = ZM0(7:11,:) + 2*kcalref*bcyc(3)*bcal(2)*0.2*abs(y(7:11,:));
ZM2(7:11,:) = ZM0(7:11,:) + 2*kcalref*bcyc(3)*bcal(2)*(0.9*abs(y(7:11,:))+2.5*x(7:11,:)-0.9*2.5);
ZM = max(ZM1,ZM2) /0.04 * Cbat;

x = x(:);
y = y(:); 
z = ZM(:);
xjcal = sdpvar(length(y),24,nBESS);
Ca = sdpvar(24,nBESS);
C = [C, xjcal >=0];
for k = 1:nBESS
    C = [C, sum(xjcal(:,1,k)) == 1, x'*xjcal(:,1,k) == SOC0(k,1)+DSOC(1,k)/2, y'*xjcal(:,1,k) == Cb0(1,k)*DSOC(1,k), z'*xjcal(:,1,k) == Ca(1,k)];
end

for i = 2:24
    y1 = (0:0.02:2.5)';
    Ty1 = Ta(i) + y1.^2*R0(i)*(9.12+3.273);
    kcycht1 = (exp(-32699/8.314*(1./Ty1-1/298.15)));
    kcyclt1 = (exp(55546/8.314*(1./Ty1-1/298.15)).*exp(2.64*(y1-2.5)/2.5));
    kch = (kcycht1*khtref+kcyclt1*kltref).* abs(y1).^0.5;
    kch  = kch.*exp(-Eacal/Rg*(1./Ty1-1/Tref));
    fun1 = @(bch,y1)bch*y1;
    bch = lsqcurvefit(fun1,0,y1,kch,-10000,10000,options);

    y2 = (-2.5:0.02:0)';
    Ty2 = Ta(i) + y2.^2*R0(i)*(9.12+3.273);
    kcycht2 = (exp(-32699/8.314*(1./Ty2-1/298.15)));
    kdis = exp(-Eacal/Rg*(1./Ty2-1/Tref)).*kcycht2*khtref.* abs(y2).^0.5;
    fun2 = @(bdis,y2)bdis*y2;
    bdis = lsqcurvefit(fun2,0,y2,kdis,-10000,10000,options);
    y3 = (-2.5:0.02:2.5)';
    funcyc = @(bcyc,y)(bcyc(1)+bcyc(2)*y).*(y<=0)+(bcyc(1)+bcyc(3)*y).*(y>0);
    bcyc = lsqcurvefit(funcyc,[0,bdis,bch],y3,[kdis; kch(2:126,1)],[0,-10000,-10000],...
        [10000,10000,10000],options);
    
    [x,y] = meshgrid(0.2:0.1:0.9,(-2.5:0.5:2.5));
    T = Ta(i) + y.^2*R0(i)*(9.12+3.273);
    kcychT = (exp(-32699/8.314*(1./T-1/298.15)));
    khT = kcychT*khtref;
    kcyclT = (exp(55546/8.314*(1./T-1/298.15)).*exp(2.64*(y-2.5)/2.5));
    klT = kcyclT*kltref;
    klT(1:(length(klT)+1)/2,:) = 0;
    kcycT = khT+klT;
    
    ZM = ((bcal(1)+bcal(2)*x(1:6,:)).*exp(-Eacal/Rg*(1./T(1:6,:)-1/Tref))*kcalref).^2 ...
    + kcycT(1:6,:).^2.*(-y(1:6,:));
    ZM0 = ZM + 2*kcalref*(bcal(1)*bcyc(2)*y(1:6,:)+bcal(2)*bcyc(1)*x(1:6,:)+bcal(1)*bcyc(1));
    ZM1 = ZM0 + 2*kcalref*(-bcyc(2))*bcal(2)*0.2*abs(y(1:6,:));
    ZM2 = ZM0 + 2*kcalref*(-bcyc(2))*bcal(2)*(0.9*abs(y(1:6,:))+2.5*x(1:6,:)-0.9*2.5);

    ZM(7:11,:) = ((bcal(1)+bcal(2)*x(7:11,:)).*exp(-Eacal/Rg*(1./T(7:11,:)-1/Tref))*kcalref).^2 ...
    + kcycT(7:11,:).^2.*y(7:11,:);
    ZM0(7:11,:) = ZM(7:11,:) + 2*kcalref*(bcal(1)*bcyc(3)*y(7:11,:)+bcal(2)*bcyc(1)*x(7:11,:)+bcal(1)*bcyc(1));
    ZM1(7:11,:) = ZM0(7:11,:) + 2*kcalref*bcyc(3)*bcal(2)*0.2*abs(y(7:11,:));
    ZM2(7:11,:) = ZM0(7:11,:) + 2*kcalref*bcyc(3)*bcal(2)*(0.9*abs(y(7:11,:))+2.5*x(7:11,:)-0.9*2.5);
    ZM = max(ZM1,ZM2) /0.04 * Cbat;
    
    x = x(:);
    y = y(:); 
    z = ZM(:);
    for k = 1:nBESS
        C = [C, sum(xjcal(:,i,k)) == 1, x'*xjcal(:,i,k) == SOC(i-1,k)+DSOC(i,k)/2, y'*xjcal(:,i,k) == Cb0(i,k)*DSOC(i,k), z'*xjcal(:,i,k) == Ca(i,k)];
    end
end

ops = sdpsettings('solver','mosek','verbose',0);
ops.mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP = 1e-12;

z = ct*P(1,:)'+sum(sum(Ca))+0.03*sum(sum(Pg));

reuslt = optimize(C,z,ops);
if reuslt.problem == 0 
    disp('success');
else
    disp('wrong');
end

I_1 = value(I);
P_1 = value(P);
Q_1 = value(Q);
V_1 = value(V);
PG = value(Pg);
QG = value(Qg);
P1 = value(P(1,:)');

SOC1 = value(SOC);
Pt1 = value(Pt);
Pac1 = value(Pac);
Pcl1 = value(Pcl);
Qt1 = value(Qt);
DSOC1 = value(DSOC);

Pt2 = zeros(24,nBESS);
for k = 1:nBESS
    Pt2(1,k) = value(U0*Cb0(1,k)*DSOC(1,k) + ...
    (0.5*ku*Cb0(1,k)+Cb0(1,k)^2*R0(1,1))*(DSOC(1,k)+ku*Cb0(1,k)/2/(0.5*ku*Cb0(1,k)+Cb0(1,k)^2*R0(1,1))*SOC0(k,1))^2+...
    (ku*Cb0(1,k))^2/4/(0.5*ku*Cb0(1,k)+Cb0(1,k)^2*R0(1,1))*(0.04-0.4*SOC0(k,1)));
    for i = 2:24
        Pt2(i,k) = value(U0*Cb0(i,k)*DSOC(i,k) + ...
        (0.5*ku*Cb0(i,k)+Cb0(i,k)^2*R0(i,1))*(DSOC(i,k)+ku*Cb0(i,k)/2/(0.5*ku*Cb0(i,k)+Cb0(i,k)^2*R0(i,1))*SOC(i-1,k))^2) + ...
        (ku*Cb0(i,k))^2/4/(0.5*ku*Cb0(i,k)+Cb0(i,k)^2*R0(i,1))*(0.04-0.4*SOC(i-1,k));
    end
end
Pcl2 = l0*Pr + l2/Pr*Pac1.^2;
Pac2 = Pcl2 + Pt2;

delta = zeros(59,24);
for j = 1:24
    for i = 1:59
        delta(i,j) = (P_1(i,j)^2+Q_1(i,j)^2)-V_1(mpc.branch(i,1),j)*I_1(i,j);
    end
end

Ca1 = value(Ca);
Ty2 = zeros(24,nBESS);
Ich2 = zeros(24,nBESS);
Cacyc2 = zeros(24,nBESS);
Cacal2 = zeros(24,nBESS);
Ca2 = zeros(24,nBESS);
Ca3 = zeros(24,nBESS);

for k = 1:nBESS
    T1 = Ta(1) + (DSOC1(1,k)*Cb0(1,k))^2*R0(1)*(9.12+3.273);
    kcalt = exp(-20592/8.314*(1/T1-1/298.15));
    xa = 8.5e-3 + (SOC0(k,1)+DSOC1(1,k)/2)*(0.78 - 8.5e-3);
    Ua = 0.6379 + 0.5416*exp(-305.5309*xa) + 0.044*tanh(-(xa-0.1958)/0.1088) ...
    - 0.1978*tanh((xa-1.0571)/0.0854) - 0.6875*tanh((xa+0.0117)/0.0529) ...
    - 0.0175*tanh((xa-0.5692)/0.0875);
    kcals = exp(0.384*96485/8.314*(0.123-Ua)/298.15) + 0.142;
    
    kcychT = (exp(-32699/8.314*(1/T1-1/298.15)));
    khT = kcychT*khtref;
    kcyclT = (exp(55546/8.314*(1/T1-1/298.15)).*exp(2.64*(DSOC1(1,k)*Cb0(1,k)-2.5)/2.5));
    klT = kcyclT*kltref;
    if DSOC1(1,k) > 0
        kcyc = khT + klT;
    else
        kcyc = khT;
    end
    Ca3(1,k) = (kcalt*kcals*kcalref + kcyc*abs(DSOC1(1,k)*Cb0(1,k))^0.5).^2 /0.04 * Cbat;
end
for j = 2:24
    for k = 1:nBESS
        T1 = Ta(j) + (DSOC1(j,k)*Cb0(j,k))^2*R0(j)*(9.12+3.273);
        kcalt = exp(-20592/8.314*(1/T1-1/298.15));
        xa = 8.5e-3 + (SOC1(j-1,k)+DSOC1(j,k)/2)*(0.78 - 8.5e-3);
        Ua = 0.6379 + 0.5416*exp(-305.5309*xa) + 0.044*tanh(-(xa-0.1958)/0.1088) ...
        - 0.1978*tanh((xa-1.0571)/0.0854) - 0.6875*tanh((xa+0.0117)/0.0529) ...
        - 0.0175*tanh((xa-0.5692)/0.0875);
        kcals = exp(0.384*96485/8.314*(0.123-Ua)/298.15) + 0.142;
    
        kcychT = (exp(-32699/8.314*(1/T1-1/298.15)));
        khT = kcychT*khtref;
        kcyclT = (exp(55546/8.314*(1/T1-1/298.15)).*exp(2.64*(DSOC1(j,k)*Cb0(j,k)-2.5)/2.5));
        klT = kcyclT*kltref;
        if DSOC1(j,k) > 0
            kcyc = khT + klT;
        else
            kcyc = khT;
        end
        Ca3(j,k) = (kcalt*kcals*kcalref + kcyc*abs(DSOC1(j,k)*Cb0(j,k))^0.5).^2 /0.04 * Cbat;
    end
end

delta1 = max(max(abs(delta)));
delta2 = max(max(abs(Pt1-Pt2)));
delta3 = max(max(abs(Pac1-l0*Pr-l2/Pr*Pac1.^2-Pt1)));


