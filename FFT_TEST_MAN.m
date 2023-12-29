function [SNR,SNDR,SFDR,THD,ENOB,FLOOR_NOISE,P_S,P_ND,fund_ind]=FFT_TEST_MAN(Vin,fs,num_H,wid,Nsample,En_plot,osr)
%功能说明
%输入输出端口定义
%输入端口
%        Vin：    需要fft的量化信号
%        fs：     adc的采样频率
%        N：      adc的位数
%        num_H:   需要计算的谐波数量
%        wid：    是否需要加窗
%        Nsample：fft采样点数
%        En_plot: 是否plot fft曲线
data=Vin;
% *********************************************
% GET DATA (These need to be set)
% ******************************
% choose the data_type which should be processed
% window or not
if  wid==0                           
    fundpnts=3;                % # of FFT points on either side of fundamental to include
                               %   in fundamental power due to non-coherence
    s_point=eps;
    windpnts=1;                % # of FFT points on either side of harmonic to include 
                           %   in harmonic power due to windowing
else
    fundpnts=3;           %   selection of fundation power
    windpnts=3;
    s_point=eps;
end
% ************************************
% SELECT WINDOW OF DATA FROM DATA SET
% ************************************
y=data(1:Nsample);
fres=fs/Nsample;
f=(0:1:Nsample-1)*fres;

% WINDOW DATA
% ***********************************
yw=y;
if wid==1;
   yw=y'.*window(@blackman,Nsample);
end
% CALCULATE VALUES FOR FFT
% ***********************************
Nsamples=round(fs/fres);   % Determine number of samples required for desired FFT
fftout = abs(fft(yw))/(Nsamples/2);            % multiply x2 to normalize because 
                                               % fft of 1*sin(x) = 0.5[d(w-wo)+d(w+wo)]
fftoutdB = 20*log10(fftout+s_point);           % change to [dB], add 1e-20 to avoid log 0


fres=fs/Nsample;  % Desired frequency resolution of FFT [Hz], fres=fclk/2^N=fin/M


fh=Nsample/osr/2;                           % frequency resolution
[Fin_dB,Fin_ind]=max(fftout(4:Nsample/2)); 
[Fin1_dB,Fin1_ind]=max(fftoutdB(4:Nsample/2)); 
fin=(Fin_ind+3-1)*fres;

% DETERMINE HARMONIC FREQUENCIES
% ***********************************
% Finds the harmoic frequencies as they fold into the baseband
%   as well as the vector index of those harmonics
for i=1:1:(num_H+1)   % harmonic number
    done=0;
    k=1;     % frequency band number, 1 = 0 to Fs Hz
    while done==0
        if i*fin < (k-1/2)*fs
            harm(i)=i*fin-(k-1)*fs;
            done=1;    
        elseif i*fin < k*fs
            harm(i)=k*fs-i*fin;
            done=1;
        else
            k=k+1;
        end
    end
end
harm_ind=harm/fres+ones(1,length(harm)); % determine vector index of harmonics
                                         % add 1 because DC = index 1
harm_ind=harm_ind';
win_ind_min=harm_ind(2:num_H+1)-windpnts;    % Points to include in harm
win_ind_max=harm_ind(2:num_H+1)+windpnts;
HD=zeros(num_H,1);

for i=1:1:num_H
    HD(i)=sum(fftout(win_ind_min(i):1:win_ind_max(i)).^2);
end

% CALCULATE SNDR, SNR, SFDR, THD
% *************************************

    % Calculate SNR/SNDR/SFDR/THD normally
    fftout_n=fftout;                                % Save FFT data, use copy for manipulations
    P_DC=sum(fftout_n(1:3).^2);
    fund_ind=(harm_ind(1)-fundpnts)...              % Points to include in fundemental
        :1:(harm_ind(1)+fundpnts);
    P_S=sum(fftout_n(fund_ind).^2);                 % Power of fundemental
 P_ND=sum(fftout_n(1:fh).^2)-P_S-P_DC;        % Power of Noise and Distortion (exclude DC)
    P_D=sum(HD(1:num_H)); % Power of harmonics
    
%    for i=[1 fund_ind]
 %       fftout_n(i)=1e-20;                            % Remove DC and fundamental from spectrum for SFDR
 %   end
    for i=1:1:3
        fftout_n(i)=1e-20;
    end
    [M_H,H_ind]=max(HD(1:num_H));        % Magnitude and index of dominant harmonic
    P_H=   M_H;                             % Power of dominant harmonic
        

    SNDRo=10*log10(P_S/P_ND);                       % SNDR [dB]
    THDo=-10*log10(P_S/P_D);                        % THD [dB]
    SNRo=10*log10(P_S/(P_ND-P_D));                  % SNR [dB]
    SFDRo=10*log10(P_S/P_H);                        % SFDR [dB]
    ENOBo=(SNDRo-1.76)/6.02;                        % ENOB [Bit]
    FLOOR_NOISEo=-SNRo-10*(log10(Nsamples/2))+fftoutdB(harm_ind(1));      % FLOOR_NOISE [dB]
    HDo=10*log10(HD(1:length(HD))/P_S);             % HD   [dB]
    SNR=SNRo;
    SFDR=SFDRo;
    
    SNDR=SNDRo;
    THD=THDo;
    HD=HDo;
    FLOOR_NOISE=FLOOR_NOISEo;
    ENOB=ENOBo;
    
% PLOT DATA, HISTOGRAM, FFT
% *********************************************
%plot FFT
%subplot(3,1,3)
if(En_plot==1)
  fffff=  fftoutdB(1:Nsamples/2)-Fin1_dB;
figure
hold on
plot(f(1:Nsamples/2)/1e6,fftoutdB(1:Nsamples/2)-Fin1_dB,'black');  % Choose FFT with frequency or index


for i=1:1:length(harm)     % Mark all the harmonics
plot(harm(i)/1e6,fftoutdB(harm_ind(i))-Fin1_dB,'r*')
text(harm(i)/1e6,fftoutdB(harm_ind(i))-Fin1_dB+3,num2str(i),'color','m','FontSize',10);
end
ylabel('Full-Scale Normalized Magnitude[dB]')
xlabel('Frequency [MHz]')
%title(sprintf('FFT (%g points)\nFs = %g MSps, Fin = %g MHz (%1.2g dBfs)\n%gx Decimated', ...
%     Nsamples,fs_predec/1e6,fin/1e6,fftoutdB(harm_ind(1)),decimate));
title(sprintf('Number of Data Points:(%g points)\n', Nsamples)); 
     
title(sprintf('Fs = %g MSps, Fin = %g MHz ', fs/1e6,fin/1e6));   %,fftoutdB(harm_ind(1))  (%1.2g dBfs)注意fin是从fft结果中计算出来的，而不是直接输入的。
      
grid;
box on;
ylim([-150 10]);
xlim([0 f(Nsamples/2)/1e6]);
set(gca,'xgrid', 'off');
set(gca, 'GridLineStyle' ,'-');
set(gca,'yTick',[-150:10:10]);

s1=sprintf('SFDR = %4.2fdB\n',SFDRo);
s2=sprintf('THD = %4.2fdB\n',THDo);
s3=sprintf('SNR   = %4.2fdB\n',SNRo);
s4=sprintf('SNDR = %4.2fdB\n',SNDRo);
s5=sprintf('ENOB = %4.2fbit\n',ENOBo);
if harm(1)/1e6 < f(Nsamples/2)/1e6/4
    xstation= f(Nsamples/2)/1e6/2;
elseif harm(1)/1e6 > (f(Nsamples/2)*3)/1e6/4
    xstation= f(Nsamples/2)/1e6/4;
else
    xstation= f(Nsamples/2)/1e6/32;
end
text(xstation,-10,s1);
text(xstation,-20,s2);
text(xstation,-30,s3);
text(xstation,-40,s4);
text(xstation,-50,s5);
hold off;  
end





if En_plot==1 % Enable plot&display function
    % % Display the Results
    % disp('DC:');
    % disp(DCsignal);
    % 
    % disp('Signal power:');
    % disp(Asignal);
    % 
    % disp('Anoise:');        
    % disp(Anoise);
    % 
    % disp('Aharmonics:');        
    % disp(Aharmo);
    % 
    % disp('HD2~HD30 Freq.(MHz):');
    % disp(FreqHD');
    % disp('HD2~HD30: HD (dBc):');
    % disp(HD');
    % % disp('Ftone:');
    % % disp(Ftone');
    % % disp('Tone:');
    % % disp(Tone');
    % 
    % disp('SFDR:');
    % disp(SFDR);
    
    disp('THD:');       
    disp(THD);

    disp('SNR:');       
    disp(SNR);

    disp('SNDR:');
    disp(SNDR);
    % 
    disp('ENOB:');
    disp(ENOB);
end
