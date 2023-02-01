function [T1,SOC1,U1,Iave,kcal,kcycht,kcyclt,kcyclths] = ...
    E_T_C_aging(U10,SOC0,T0,Tf,Pt,Cb)

al = 1/(9.12+3.273);
C = 73.2;
D = 3600*Cb;


opts = odeset('RelTol',1e-2,'AbsTol',1e-5,'MaxStep',1);
[t,y] = ode45(@(t,y) odefun(t,y,al,C,D,Pt,Tf),[0 3600],[T0; SOC0; U10],opts);
% [t,y] = ode45(@(t,y) odefun(t,y,al,C,D,Pt,Tf),[0 3600],[T0; SOC0; U10]);
T1 = y(length(y(:,1)),1);
SOC1 = y(length(y(:,1)),2);
U1 = y(length(y(:,1)),3);

I = (((Uocv(y(:,2))+y(:,3)).^2+4*Pt.*(0.0976*exp(-0.0077*y(:,1)))).^0.5-(Uocv(y(:,2))+y(:,3)))/2./(0.0976*exp(-0.0077*y(:,1)));
Iave = mean(I);

Tc = y(:,1) + 273.15;
kltref = 4.01e-4;
khtref = 1.46e-4;
klthsref = 2.03e-6;
kcalref = 3.69e-4 * 1.5;

a = 3;
kltref = a*kltref;
khtref = a*khtref;
klthsref = a*klthsref;

% 日历老化
kcalt = exp(-20592/8.314*(1./Tc-1/298.15));
xa = 8.5e-3 + y(:,2)*(0.78 - 8.5e-3);
Ua = 0.6379 + 0.5416*exp(-305.5309*xa) + 0.044*tanh(-(xa-0.1958)/0.1088) ...
    - 0.1978*tanh((xa-1.0571)/0.0854) - 0.6875*tanh((xa+0.0117)/0.0529) ...
    - 0.0175*tanh((xa-0.5692)/0.0875);
kcals = exp(0.384*96485/8.314*(0.123-Ua)/298.15) + 0.142;
kcal = mean(kcalt.*kcals) * kcalref;

% 循环老化
% I = I/2.5*3;
Ich = max([I zeros(length(I),1)]')';
klt = exp(55546/8.314*(1./Tc-1/298.15)).*exp(2.64*(Ich-2.5)/2.5);
klt(Ich == 0,1) = 0;
kht = exp(-32699/8.314*(1./Tc-1/298.15));
klths = exp(2.3e5/8.314*(1./Tc-1/298.15)).*exp(7.8*(I-2.5)/2.5);
klths(Ich == 0,1) = 0;
kcycht = mean(kht) * khtref;
kcyclt = mean(klt) * kltref;
kcyclths = mean(klths) * klthsref;


function dydt = odefun(t,y,al,C,D,Pt,Tf)
% y1 = T; y2 = SOC; y3 = U1;
dydt = [((Uocv(y(2))+y(3))*((Uocv(y(2))+y(3))-((Uocv(y(2))+y(3))^2+4*Pt*(0.0976*exp(-0.0077*y(1))))^0.5)/2/(0.0976*exp(-0.0077*y(1)))+Pt-al*(y(1)-Tf))/C; ...
        (((Uocv(y(2))+y(3))^2+4*Pt*(0.0976*exp(-0.0077*y(1))))^0.5-(Uocv(y(2))+y(3)))/2/(0.0976*exp(-0.0077*y(1)))/D/(-0.354*exp(-0.031*y(1))+1.157); ...
        -y(3)/R1(y(1))/C1(y(1))+(-(Uocv(y(2))+y(3))+(((Uocv(y(2))+y(3)))^2+4*Pt*(0.0976*exp(-0.0077*y(1))))^0.5)/2/(0.0976*exp(-0.0077*y(1)))/C1(y(1))];
end

function U = Uocv(SOC)
    U = 3.44-0.191*SOC+0.108*log(SOC)-0.033*log(1-SOC); 
end

function R1 = R1(Tf)
    R1 = 2.5e-8*Tf^2+2e-6*Tf+6e-5; 
end

function C1 = C1(Tf)
    C1 = -0.465*Tf^3 + 42*Tf^2 -1010*Tf +9400; 
end

end





