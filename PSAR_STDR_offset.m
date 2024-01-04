clear
% 定义ADC的基本参数
Vref = 1;                                       % 定义基准电压为1V 
Vcm = 0.5;                                      % 定义共模电平为0.5V
N  = 12;                                          %总位数
N1 = 5;                                          % 第一级位数
N2 = N - N1 ;                                    % 第二级位数 

% 定义失调量
Num_os = 61;                                      % 失调仿真个数
Comp_os = linspace(0,0,Num_os);                   % 比较器失调
Amp_os =  linspace(-1e-3,1e-3,Num_os);          % 运放失调

% 定义电容失配的基本参数
sig_c1 = 0;                                      % 定义单位电容的标准偏差 
C_nor1 = 2.^[N1-1:-1:0];                            % 定义理想电容阵列, 其中默认电位电容为1 
Wda1 = 2.^[N1-1:-1:0];

sig_c2 = 0;
C_nor2 = 2.^[N2-2:-1:0];                        % 定义理想电容阵列, 其中默认电位电容为1 
Wda2 = 2.^[N2-2:-1:0];

% 产生斜坡信号，针对其中均匀分布的2N_in个电压值进行处理 
M = 5;                                            %定义每个转换台阶包含的点数为2^M 
N_in = N1 + N2 + M; 
Vip_ramp = [0:1:2^N_in-1]/2^N_in*Vref;                 % 产生斜坡输入 
Vin_ramp = Vref - Vip_ramp;

% 产生正弦信号
fclk = 100;                       % 100MHz 采样率
M_sin    = 9  ;                        % 采样时间内包含的输入信号的周期数
N_sin    = 1024;                       % 采样点数
fin=M_sin/N_sin*fclk;                     % 输入信号频率
T1=1:N_sin;
Vip_sin = 0.5 + 0.5*sin(2*pi*fin/fclk*T1);
Vin_sin = 0.5 - 0.5*sin(2*pi*fin/fclk*T1);


for os_num = 1:Num_os            % 比较器失调的个数
    % 电容失配的蒙特卡洛次数
    Num = 1;                                       
    
    % 开始转换
    for j = 1:Num
        C_dev1(j,:) = sig_c1*sqrt(C_nor1).*randn(1,N1);          % 电容阵列中各电容的标准偏差，呈正态分布
        C_act1(j,:) = C_nor1 + C_dev1(j,:);                            % 定义实际电容由其均值和标准偏差组成
        C_dev2(j,:) = sig_c2*sqrt(C_nor2).*randn(1,N2-1);          % 电容阵列中各电容的标准偏差，呈正态分布
        C_act2(j,:) = C_nor2 + C_dev2(j,:);                            % 定义实际电容由其均值和标准偏差组成
        for i = 1:2^N_in 
            % 量化斜坡信号
            [D1_ramp(j,i),Vres_p_ramp(j,i),Vres_n_ramp(j,i)] = Coarse_sar(Vip_ramp(i), Vin_ramp(i), Vref, Vcm, N1, C_act1(j,:), C_act1(j,:) , 1, 1, 0, 0, Comp_os(os_num), 0, 0, Wda1);                  % 针对均匀分布的2N_in个电压值进行A/D转换 

            V_residue_ramp(j,i) = (Vres_n_ramp(j,i) - Vres_p_ramp(j,i)+ Amp_os(os_num) )*2^N1;
            V_residue_p_ramp(j,i) = Vcm + V_residue_ramp(j,i)/2;
            V_residue_n_ramp(j,i) = Vcm - V_residue_ramp(j,i)/2;

            D2_ramp(j,i) = Fine_sar(V_residue_p_ramp(j,i), V_residue_n_ramp(j,i), Vref, Vcm, N2-1, C_act2(j,:), C_act2(j,:) , 1, 1, 0, 0, 0, 0, 0, Wda2);
            D_ramp(j,i) = D1_ramp(j,i)*2^N2 + D2_ramp(j,i);
            Do_ramp(j,i) = D_ramp(j,i)*2*Vref/(2^(N1+N2))-1;      % 输出差分信号范围[-Vref,Vref]
        end

        % 画余差曲线
         hold on
        plot(Vip_ramp, V_residue_p_ramp);
        xlabel('V_{idiff}'); % 添加横坐标标识
        ylabel('V_{res}'); % 添加纵坐标标识
        title('Residue Characteristic Curve'); % 添加标题
         hold off

        for i = 1:N_sin
            % 量化正弦信号
            [D1_sin(j,i),Vres_p_sin(j,i),Vres_n_sin(j,i)] = Coarse_sar(Vip_sin(i), Vin_sin(i), Vref, Vcm, N1, C_act1(j,:), C_act1(j,:) , 1, 1, 0, 0, Comp_os(os_num), 0, 0, Wda1);                  % 针对均匀分布的2N_in个电压值进行A/D转换 

            V_residue_sin(j,i) = (Vres_n_sin(j,i) - Vres_p_sin(j,i)+ Amp_os(os_num) )*2^N1;
            V_residue_p_sin(j,i) = Vcm + V_residue_sin(j,i)/2;
            V_residue_n_sin(j,i) = Vcm - V_residue_sin(j,i)/2;

            D2_sin(j,i) = Fine_sar(V_residue_p_sin(j,i), V_residue_n_sin(j,i), Vref, Vcm, N2-1, C_act2(j,:), C_act2(j,:) , 1, 1, 0, 0, 0, 0, 0, Wda2);
            D_sin(j,i) = D1_sin(j,i)*2^N2 + D2_sin(j,i);
            Do_sin(j,i) = D_sin(j,i)*2*Vref/(2^(N1+N2))-1;      % 输出差分信号范围[-Vref,Vref]
        end 
    end






    % 计算DNL的标准差
    for j = 1:Num
        for i = 1:2^(N1+N2) 
            A = (i-1)*ones(1,2^N_in); 
            E = D_ramp(j,:) - A; 
            Q(i) = 2^N_in - nnz(E); 
            DNL(j,i)=Q(i)/2^M - 1; 
            INL(j,i) = sum(DNL(j,(1:i)));   % 计算A/D转换器的DNL和INL 
        end 

    end

    DNL_mean = mean(DNL,2);                  %正常情况下均值应该是0

    DNL_std = std(DNL,0,2);

    DNL_final(os_num) = mean(DNL_std);

    % FFT 计算ENOB
    num_H=5;
    wid=0;
    En_plot=0;
    osr=1;
    fs=100*10^6;
    Nsample=1024;
    Do_sin1 = Do_sin';
    for j = 1:Num
    [SNR(j),SNDR(j),SFDR(j),THD(j),ENOB(j),FLOOR_NOISE(j),P_S(j),P_ND(j),fund_ind]=FFT_TEST_MAN(Do_sin1(:,j),fs,num_H,wid,Nsample,En_plot,osr);
    end

    ENOB_final(os_num) = mean(ENOB);
end
