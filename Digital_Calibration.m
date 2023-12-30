clear
Bc = zeros(1,4);
B = zeros(1, 9);
B1 = [B,Bc];

C = [254.38,122,68,36,22,14,8,4,2,1]; %最高可对278的权重校准
Cc = [4,2,1,1];
Wc = [1/2,1/4,1/8,1/8];

Cu = 1;
Ccu = 1;
Vcm = 0.5;
vref = 1;
csum = sum(C) + Cu + Ccu;

Vcalp = Vcm + B1(1:4) * Wc' * vref ;
Vcaln = Vcm + B1(1:4) * Wc' * vref ;

Vdacp = Vcm - ((C(1) - B1(1:9)*C(2:10)')*vref + (Vcalp - Vcm)/Ccu)/csum;
Vdacn = Vcm + ((C(1) - B1(1:9)*C(2:10)')*vref + (Vcalp - Vcm)/Ccu)/csum;

for i = 1:9
    B1(i) = 1;
    Vdacp = Vcm - ((C(1) - B1(1:9)*C(2:10)')*vref + (Vcalp - Vcm)/Ccu)/csum;
    Vdacn = Vcm + ((C(1) - B1(1:9)*C(2:10)')*vref + (Vcalp - Vcm)/Ccu)/csum;
    if (Vdacp <= Vdacn)
        B1(i) = 1;
    else
        B1(i) = 0;
    end
end

Vdacp = Vcm - ((C(1) - B1(1:9)*C(2:10)')*vref - (Vcalp - Vcm)/Ccu)/csum;
Vdacn = Vcm + ((C(1) - B1(1:9)*C(2:10)')*vref - (Vcalp - Vcm)/Ccu)/csum;

for i = 10 : 13
    B1(i) = 1;
    Vcalp = Vcm + B1(10:13) * Wc' * vref ;
    Vcaln = Vcm + B1(10:13) * Wc' * vref ;
    
    Vdacp = Vcm - ((C(1) - B1(1:9)*C(2:10)')*vref - (Vcalp - Vcm)/Ccu)/csum;
    Vdacn = Vcm + ((C(1) - B1(1:9)*C(2:10)')*vref - (Vcalp - Vcm)/Ccu)/csum;
    if (Vdacp <= Vdacn)
        B1(i) = 1;
    else
        B1(i) = 0;
    end
    
end


w = B1(1:9) * C(2:10)' + B1(10:13) * Wc'