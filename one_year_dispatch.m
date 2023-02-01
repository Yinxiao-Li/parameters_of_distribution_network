load weather;
load kload;
SB = 1;
Cba = 2.5*ones(1,5);
Cbn = 2.5;
Cbat = 1800;

SOC0a = 0.2*ones(5,1);
T0a = Tf(1,1)*ones(5,1);
U1a0 = zeros(5,1);

Callossa = zeros(24,5);
Cyclosshta = zeros(24,5);
Cyclosslta = zeros(24,5);
Cyclosslthsa = zeros(24,5);
tacc = 0;
Qacc = zeros(1,5);
Qloss = zeros(366,5);

SOC8760 = zeros(8760,5);
Tb8760 = zeros(8760,5);
U18760 = zeros(8760,5);
kcala8760 = zeros(8760,5);
kcychta8760 = zeros(8760,5);
kcyclta8760 = zeros(8760,5);
kcyclthsa8760 = zeros(8760,5);
Iavea8760 = zeros(8760,5);
Qacc8760 = zeros(8761,5);

delta365 = zeros(365,3);

Cb365 = zeros(365,5);
Cgridact365 = zeros(365,1); 
Cgridexp365 = zeros(365,1); 
Caact365 = zeros(365,1); 
Caexp365 = zeros(365,1); 

PBEact365 = zeros(365,1); 
PBEexp365 = zeros(365,1); 

for j = 1:365
    tic
    if j > 151 && j < 244
        ct = 1000*SB*[0.353 0.353 0.353 0.353 0.353 0.353 0.742 0.742 1.196 1.196 1.196 0.742 0.742 1.196 1.196 0.742 0.742 0.742 1.196 1.196 1.196 0.742 0.353 0.353];
    else
        ct = 1000*SB*[0.353 0.353 0.353 0.353 0.353 0.353 0.742 0.742 1.196 1.196 1.196 0.742 0.742 0.742 0.742 0.742 0.742 0.742 1.196 1.196 1.196 0.742 0.353 0.353];
    end
    kload24 = kload(24*j-23:24*j,1);
    Tf24 = Tf(24*j-23:24*j,1);
    W24 = W(24*j-23:24*j,1);
    
    SOC0a = 0.2*ones(5,1);
    
    % solve the dispatch model
    [P1a, Pt1a, Ca1, PGa, QGa, Pac1a, Qt1a, SOC1a, delta1, delta2, delta3] = opt_oneday(SOC0a, W24, Tf24, ct, kload24, Cba);
    delta365(j,1) = delta1;
    delta365(j,2) = delta2;
    delta365(j,3) = delta3;
    
    % calculate the actual power of BESSs
    Tb_a = zeros(24,5);
    SOCa = zeros(24,5);
    U1a = zeros(24,5);
    Cratea = zeros(24,5);
    Pt15a = Pt1a;

    for k = 1:5
        [Tb_a(1,k),SOCa(1,k),U1a(1,k),Cratea(1,k)] = E_T_C(U1a0(k,1),SOC0a(k,1),T0a(k,1),Tf24(1),Pt15a(1,k),Cba(1,k));
        v1 = 0;
        v2 = 0;
        while Cratea(1,k) > 0.5 || SOCa(1,k)>0.9 || SOCa(1,k)<0.2
            Pt15a(1,k) = Pt15a(1,k) - 0.1;
            [Tb_a(1,k),SOCa(1,k),U1a(1,k),Cratea(1,k)] = E_T_C(U1a0(k,1),SOC0a(k,1),T0a(k,1),Tf24(1),Pt15a(1,k),Cba(1,k));
            v1 = 1;
        end
        if v1 > 0
            Pt15a(1,k) = Pt15a(1,k) + 0.1;
            [Tb_a(1,k),SOCa(1,k),U1a(1,k),Cratea(1,k)] = E_T_C(U1a0(k,1),SOC0a(k,1),T0a(k,1),Tf24(1),Pt15a(1,k),Cba(1,k));
        end
        while Cratea(1,k) > 0.5 || SOCa(1,k)>0.9 || SOCa(1,k)<0.2
            Pt15a(1,k) = Pt15a(1,k) - 0.01;
            [Tb_a(1,k),SOCa(1,k),U1a(1,k),Cratea(1,k)] = E_T_C(U1a0(k,1),SOC0a(k,1),T0a(k,1),Tf24(1),Pt15a(1,k),Cba(1,k));
            v2 = 1;
        end
        if v2 > 0
            Pt15a(1,k) = Pt15a(1,k) + 0.01;
            [Tb_a(1,k),SOCa(1,k),U1a(1,k),Cratea(1,k)] = E_T_C(U1a0(k,1),SOC0a(k,1),T0a(k,1),Tf24(1),Pt15a(1,k),Cba(1,k));
        end
        while Cratea(1,k) > 0.5 || SOCa(1,k)>0.9 || SOCa(1,k)<0.2
            Pt15a(1,k) = Pt15a(1,k) - 0.001;
            [Tb_a(1,k),SOCa(1,k),U1a(1,k),Cratea(1,k)] = E_T_C(U1a0(k,1),SOC0a(k,1),T0a(k,1),Tf24(1),Pt15a(1,k),Cba(1,k));
        end
    end

    for k = 1:5
        for i = 2:24
            [Tb_a(i,k),SOCa(i,k),U1a(i,k),Cratea(i,k)] = E_T_C(U1a(i-1,k),SOCa(i-1,k),Tb_a(i-1,k),Tf24(i),Pt15a(i,k),Cba(1,k));
            v1 = 0;
            v2 = 0;
            Pt00 = Pt15a(i,k);
            while Cratea(i,k) > 0.5 || SOCa(i,k)>0.9 || SOCa(i,k)<0.2
                if Pt00 > 0
                    Pt15a(i,k) = Pt15a(i,k) - 0.1;
                else
                    Pt15a(i,k) = Pt15a(i,k) + 0.1;
                end
                [Tb_a(i,k),SOCa(i,k),U1a(i,k),Cratea(i,k)] = E_T_C(U1a(i-1,k),SOCa(i-1,k),Tb_a(i-1,k),Tf24(i),Pt15a(i,k),Cba(1,k));
                v1 = 1;
            end
            if v1 > 0
                if Pt00 > 0
                    Pt15a(i,k) = Pt15a(i,k) + 0.1;
                else
                    Pt15a(i,k) = Pt15a(i,k) - 0.1;
                end
                [Tb_a(i,k),SOCa(i,k),U1a(i,k),Cratea(i,k)] = E_T_C(U1a(i-1,k),SOCa(i-1,k),Tb_a(i-1,k),Tf24(i),Pt15a(i,k),Cba(1,k));
            end
            while Cratea(i,k) > 0.5 || SOCa(i,k)>0.9 || SOCa(i,k)<0.2
                if Pt15a(i,k) > 0
                    Pt15a(i,k) = Pt15a(i,k) - 0.01;
                else
                Pt15a(i,k) = Pt15a(i,k) + 0.01;
                end
                v2 = 1;
                [Tb_a(i,k),SOCa(i,k),U1a(i,k),Cratea(i,k)] = E_T_C(U1a(i-1,k),SOCa(i-1,k),Tb_a(i-1,k),Tf24(i),Pt15a(i,k),Cba(1,k));
            end
            if v2 > 0
                if Pt00 > 0
                    Pt15a(i,k) = Pt15a(i,k) + 0.01;
                else
                    Pt15a(i,k) = Pt15a(i,k) - 0.01;
                end
                [Tb_a(i,k),SOCa(i,k),U1a(i,k),Cratea(i,k)] = E_T_C(U1a(i-1,k),SOCa(i-1,k),Tb_a(i-1,k),Tf24(i),Pt15a(i,k),Cba(1,k));
            end
            while Cratea(i,k) > 0.5 || SOCa(i,k)>0.9 || SOCa(i,k)<0.2
                if Pt15a(i,k) > 0
                    Pt15a(i,k) = Pt15a(i,k) - 0.001;
                else
                    Pt15a(i,k) = Pt15a(i,k) + 0.001;
                end
                [Tb_a(i,k),SOCa(i,k),U1a(i,k),Cratea(i,k)] = E_T_C(U1a(i-1,k),SOCa(i-1,k),Tb_a(i-1,k),Tf24(i),Pt15a(i,k),Cba(1,k));
            end
        end
    end
    
    for k = 1:5
        for i = 24:24
            [Tb_a(i,k),SOCa(i,k),U1a(i,k),Cratea(i,k)] = E_T_C(U1a(i-1,k),SOCa(i-1,k),Tb_a(i-1,k),Tf24(i),Pt15a(i,k),Cba(1,k));
            v1 = 0;
            v2 = 0;
            Pt00 = Pt15a(i,k);
            while SOCa(i,k)>0.2
                Pt15a(i,k) = Pt15a(i,k) - 0.1;
                [Tb_a(i,k),SOCa(i,k),U1a(i,k),Cratea(i,k)] = E_T_C(U1a(i-1,k),SOCa(i-1,k),Tb_a(i-1,k),Tf24(i),Pt15a(i,k),Cba(1,k));
            end
            Pt15a(i,k) = Pt15a(i,k) + 0.1;
            [Tb_a(i,k),SOCa(i,k),U1a(i,k),Cratea(i,k)] = E_T_C(U1a(i-1,k),SOCa(i-1,k),Tb_a(i-1,k),Tf24(i),Pt15a(i,k),Cba(1,k));
            
            while SOCa(i,k)>0.2
                Pt15a(i,k) = Pt15a(i,k) - 0.01;
                [Tb_a(i,k),SOCa(i,k),U1a(i,k),Cratea(i,k)] = E_T_C(U1a(i-1,k),SOCa(i-1,k),Tb_a(i-1,k),Tf24(i),Pt15a(i,k),Cba(1,k));
            end
            Pt15a(i,k) = Pt15a(i,k) + 0.01;
            [Tb_a(i,k),SOCa(i,k),U1a(i,k),Cratea(i,k)] = E_T_C(U1a(i-1,k),SOCa(i-1,k),Tb_a(i-1,k),Tf24(i),Pt15a(i,k),Cba(1,k));
            
            while SOCa(i,k)>0.2
                Pt15a(i,k) = Pt15a(i,k) - 0.001;
                [Tb_a(i,k),SOCa(i,k),U1a(i,k),Cratea(i,k)] = E_T_C(U1a(i-1,k),SOCa(i-1,k),Tb_a(i-1,k),Tf24(i),Pt15a(i,k),Cba(1,k));
            end
            Pt15a(i,k) = Pt15a(i,k) + 0.001;
            [Tb_a(i,k),SOCa(i,k),U1a(i,k),Cratea(i,k)] = E_T_C(U1a(i-1,k),SOCa(i-1,k),Tb_a(i-1,k),Tf24(i),Pt15a(i,k),Cba(1,k));
            
        end
    end
    
    % calculate the actual aging of BESSs
    Tb_a = zeros(24,5);
    SOCa = zeros(24,5);
    U1a = zeros(24,5);
    Iavea = zeros(24,5);
    kcala = zeros(24,5);
    kcychta = zeros(24,5);
    kcyclta = zeros(24,5);
    kcyclthsa = zeros(24,5);
    for k = 1:5
        [Tb_a(1,k),SOCa(1,k),U1a(1,k),Iavea(1,k),kcala(1,k),kcychta(1,k),kcyclta(1,k),...
        kcyclthsa(1,k)] = E_T_C_aging(U1a0(k,1),SOC0a(k,1),T0a(k,1),Tf24(1),Pt15a(1,k),Cba(1,k));
        for i = 2:24
            [Tb_a(i,k),SOCa(i,k),U1a(i,k),Iavea(i,k),kcala(i,k),kcychta(i,k),kcyclta(i,k),kcyclthsa(i,k)]...
            = E_T_C_aging(U1a(i-1,k),SOCa(i-1,k),Tb_a(i-1,k),Tf24(i),Pt15a(i,k),Cba(1,k));
        end
    end

    for k = 1:5
        for i = 1:24
            Callossa(i,k) = kcala(i,k) * ((tacc+1)^0.5-tacc^0.5);
            Cyclosshta(i,k) = kcychta(i,k) * ((Qacc(1,k)+abs(Iavea(i,k)))^0.5-Qacc(1,k)^0.5);
            Cyclosslta(i,k) = kcyclta(i,k) * ((Qacc(1,k)+abs(Iavea(i,k)))^0.5-Qacc(1,k)^0.5);
            Cyclosslthsa(i,k) = kcyclthsa(i,k) * max([0 SOCa(i,k)-0.82]) * Iavea(i,k);
            tacc = tacc + 1;
            Qacc(1,k) = Qacc(1,k) + abs(Iavea(i,k));
            Qacc8760(24*j-23+i,k) = Qacc8760(24*j-24+i,k) + abs(Iavea(i,k));
        end
        tacc = tacc - 24;
    end
    tacc = tacc + 24;
    Qloss(j+1,:) = Qloss(j,:) + sum(Callossa)+sum(Cyclosshta)+sum(Cyclosslta)+sum(Cyclosslthsa);
    % actual aging cost
    Caging2a = zeros(1,5);
    for k = 1:5
        Caging2a(1,k) = sum((kcala(:,k) + (kcychta(:,k)+kcyclta(:,k)).*abs(Iavea(:,k)).^0.5).^2/0.04*2.5*3.3*2424/1000*Cbat);
    end

    % calculate the actual purchase cost
    l0 = 7.503001200480192e-04;
    l2 = 0.025510204081633;
    Pr = 2.5*3.3; 

    Pac15a = (1 - (1 - 4*l2/Pr*(l0*Pr + Pt15a)).^0.5)/2/(l2/Pr);

    load distribution_network;
    SB = mpc.baseMVA;
    Pz = mpc.bus(1:60,3);
    Qz = mpc.bus(1:60,4);
    Pz = Pz / SB;
    Qz = Qz / SB;

    nbat = 2424;
    P11 = zeros(24,1);
    for i = 1:24
        mpc.bus(:,3) = kload24(i)*Pz * SB;
        mpc.bus(:,4) = kload24(i)*Qz * SB;
        mpc.bus(42,3) = kload24(i)*Pz(42) * SB - PGa(1,i) * SB + Pac15a(i,1)*nbat/1000000;
        mpc.bus(42,4) = kload24(i)*Qz(42) * SB - QGa(1,i) * SB + Qt1a(i,1)*nbat/1000000;
        mpc.bus(48,3) = kload24(i)*Pz(48) * SB - PGa(2,i) * SB + Pac15a(i,2)*nbat/1000000;
        mpc.bus(48,4) = kload24(i)*Qz(48) * SB - QGa(2,i) * SB + Qt1a(i,2)*nbat/1000000;
        mpc.bus(32,3) = kload24(i)*Pz(32) * SB - PGa(3,i) * SB + Pac15a(i,3)*nbat/1000000;
        mpc.bus(32,4) = kload24(i)*Qz(32) * SB - QGa(3,i) * SB + Qt1a(i,3)*nbat/1000000;
        mpc.bus(45,3) = kload24(i)*Pz(45) * SB - PGa(4,i) * SB + Pac15a(i,4)*nbat/1000000;
        mpc.bus(45,4) = kload24(i)*Qz(45) * SB - QGa(4,i) * SB + Qt1a(i,4)*nbat/1000000;
        mpc.bus(47,3) = kload24(i)*Pz(47) * SB - PGa(5,i) * SB + Pac15a(i,5)*nbat/1000000;
        mpc.bus(47,4) = kload24(i)*Qz(47) * SB - QGa(5,i) * SB + Qt1a(i,5)*nbat/1000000;
        mpopt = mpoption('out.all', 0 ,'verbose', 0);
        result1 = runpf(mpc,mpopt);
        P11(i,1) = result1.branch(1,14);
    end

    cta = ct * P11;
    
    SOC8760(24*j-23:24*j,:) = SOCa;
    Tb8760(24*j-23:24*j,:) = Tb_a;
    U18760(24*j-23:24*j,:) = U1a;
    kcala8760(24*j-23:24*j,:) = kcala;
    kcychta8760(24*j-23:24*j,:) = kcychta;
    kcyclta8760(24*j-23:24*j,:) = kcyclta;
    kcyclthsa8760(24*j-23:24*j,:) = kcyclthsa;
    Iavea8760(24*j-23:24*j,:) = Iavea;
    
    Cb365(j,:) = Cbn*(1-Qloss(j+1,:));
    
    SOC0a = SOCa(24,:)';
    T0a = Tb_a(24,:)';
    U1a0 = U1a(24,:)';
    Cba = Cb365(j,:);
    
    Cgridact365(j,1) = cta+30*sum(sum(PGa));
    Cgridexp365(j,1) = ct*P1a+30*sum(sum(PGa));
    Caact365(j,1) = sum(Caging2a);
    Caexp365(j,1) = sum(sum(Ca1));
    
    PBEact365(j,1) = sum(sum(abs(Pac15a))); 
    PBEexp365(j,1) = sum(sum(abs(Pac1a)));

save('result.mat')
disp(j);

toc
end
