function [T,SOC,U1,Cmax]=E_T_C(U10,SOC0,T0,Tf,Pt,Cb)

al = 1/(9.12+3.273);
C = 73.2;
D = 3600*Cb;


opts = odeset('RelTol',1e-2,'AbsTol',1e-5,'MaxStep',2);
[t,y] = ode45(@(t,y) odefun(t,y,al,C,D,Pt,Tf),[0 3600],[T0; SOC0; U10],opts);
% [t,y] = ode45(@(t,y) odefun(t,y,al,C,D,Pt,Tf),[0 3600],[T0; SOC0; U10]);
T = y(length(y(:,1)),1);
SOC = y(length(y(:,1)),2);
U1 = y(length(y(:,1)),3);
I = (((Uocv(y(:,2))+y(:,3)).^2+4*Pt.*(0.0976*exp(-0.0077*y(:,1)))).^0.5-(Uocv(y(:,2))+y(:,3)))/2./(0.0976*exp(-0.0077*y(:,1)));
Crate = I/Cb./(-0.354*exp(-0.031*y(:,1))+1.157);
Cmax = max(abs(Crate));


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





