function [adout] = Fine_sar(Vip, Vin, Vref, Vcm, N, Cm_p, Cm_n, Cd1_p, Cd1_n, Cp1_p, Cp1_n, Comp_os, del_Compvn, del_ktc, Wda)
Ctm_p = sum(Cm_p) + Cd1_p + Cp1_p;  % Cm_p 权重电容 Cd1_p单位电容 Cp1_p寄生电容
Ctm_n = sum(Cm_n) + Cd1_n + Cp1_n;
Qi_p = -(sum(Cm_p) + Cd1_p) * (Vip-Vcm+del_ktc*randn(1,1));
Qi_n = -(sum(Cm_n) + Cd1_n) * (Vin-Vcm+del_ktc*randn(1,1));
Vdam_p(1:N) = Vcm; % define bottom plate voltages
Vdam_n(1:N) = Vcm;
adout = 0; % output

for i = 1 : N
    Vres_p = Vcm + (Qi_p + (Vdam_p-Vcm)*rot90(Cm_p,3))/Ctm_p;
    Vres_n = Vcm + (Qi_n + (Vdam_n-Vcm)*rot90(Cm_n,3))/Ctm_n;
    if Vres_p-Vres_n >= Comp_os+del_Compvn*randn(1,1)
        Vdam_p(i) = 0; Vdam_n(i) = Vref;
    else
        Vdam_p(i) = Vref; Vdam_n(i) = 0;
		adout = adout + Wda(i)*2;
    end
    
end

Vres_p = Vcm + (Qi_p + (Vdam_p-Vcm)*rot90(Cm_p,3))/Ctm_p;
Vres_n = Vcm + (Qi_n + (Vdam_n-Vcm)*rot90(Cm_n,3))/Ctm_n;
if Vres_p-Vres_n <= Comp_os+del_Compvn*randn(1,1) adout = adout + 1; end