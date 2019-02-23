%2000 us de periodo de muestreo
%115200 puerto
%19,26,27,02,33,39(incompleto), NO
%38 en adelante las ultimas

clear,clc

datalb= dlmread('steering.log');
data=dlmread('steering2.log');
%24,25 MENOS DE 45000
%slb=datalb(end-45000:end,1).*5/1024;%Linea base conductividad
%plb=datalb(end-45000:end,2).*5/1024;%Linea base pulso

% slb=datalb(:,1).*5/1024;%Linea base conductividad
% plb=datalb(:,2).*5/1024;%Linea base pulso
 
%Sad
s=data(:,1).*5/1024; %Conductividad
p=data(:,2).*5/1024; %Pulso

Fs=500; %F. sampling
T=1/Fs; %T. sampling

%%
Tlb=length(slb);    %Numero de muestras de linea base
Tamp=length(p);     %Numero de muestras de tratamiento

%Vector de tiempo  tratamiento
t=(0:Tamp-1)*T; 
t1=t(1001:end); %Vector de tiempo para conductividad
t2=(0:Tamp-1)*T;%Vector de tiempo para la señal de Temperatura
t=t(51:end);    %Vector de tiempo para ritmo

%Vector de tiempo  linea base
Tbase=(0:Tlb-1)*T;
tglb=Tbase(1001:end);%Vector de tiempo para conductividad
tplb=Tbase(51:end); %Vector de tiempo para ritmo

%Definicion ventana hamming 
 
win=hamming(101);   %Tamaño de ventana para Ritmo: 101
win2=hamming(1001); %Tamaño de ventana para Conductividad: 1001

%Filtro pasabanda tipo FIR - Ritmo cardiaco
f1=0.5;
f2=10;

w1=2*f1/Fs;
w2=2*f2/Fs;
B=fir1(100,[w1 w2],'bandpass',win); %Filtro pasabanda tipo FIR
y=filter(B,1,p);
y=y(51:end);        %Eliminacion de primeras 50 muestras
pulsoaf=y;   

%Filtro pasabanda tipo FIR - Ritmo cardiaco linea base
ylb=filter(B,1,plb);
ylb=ylb(51:end);        %Eliminacion de primeras 50 muestras
pulsolb=ylb;   


%Filtrado pasabajos GSR
f3=1;
w3=2*f3/Fs;
B2=fir1(1000,w3,'low',win2);
y1=filter(B2,1,s);
y1=y1(1001:end);      %Eliminacion de primeras 1000 muestras
gsr = y1;
%Filtro de media movil
gsrmm=movmean(gsr,500); %Mejor resultado tomando el valor medio cada 500 muestras (1 segundo)


%Filtrado linea base GSR
y2=filter(B2,1,slb);
y2=y2(1001:end);      %Eliminacion de primeras 1000 muestras
gsrlb = y2;
%Filtro de media movil
gsrmmlb=movmean(gsrlb,500); %Mejor resultado tomando el valor medio cada 500 muestras (1 segundo)


%Segmentacion en ventanas de 5 segundos
 
 bx = 2500; %Division del total de los clips en ventanas de 5 segundos  
 na = numel(pulsoaf); %Numero de elementos de la señal de Ritmo cardiaco
 c = mat2cell(pulsoaf,diff([0:bx:na-1,na]));%Division de señal en ventanas de 5000 muestras
 c=c(1:end-1);
 n_iter = length(c);
 
 
%Segmentacion en ventanas de 5 segundos (Linea base)

 blb= 2500; %Division del total de los clips en ventanas de 5 segundos  
 nlb = numel(pulsolb); %Numero de elementos de la señal de Ritmo cardiaco
 clb = mat2cell(pulsolb,diff([0:blb:nlb-1,nlb]));%Division de señal en ventanas de 5000 muestras
 clb=clb(1:end-1);
 n_iterlb = length(clb);
 
%% 
%Definicion ventana hamming para FFT
winh=hamming(2500);
%Calculo de la FFT para la linea base del ritmo cardiaco
 for i=1:n_iterlb
 c_lb=clb{i}-mean(clb{i});
 Xlb=fft(clb{i}.*winh,15000);
 bpmlb=abs(Xlb);
 bpmlb(1:30)=0;                 
 bpmlb(54:end)=0;
 [ylb,in]=max(bpmlb);
 hrlb(i)=(in-1)*60*Fs/(length(Xlb)); % Vector de ritmo cardiaco de linea base en BPM
 Perlb = round(15000/(in-1));
 %extracci�n de arm�nicos
 cSFlb=abs(fft(c_lb(1:Perlb)))/ylb;
 ARMlb(i,:)=cSFlb(3:6);
 end
 
 for i=1:n_iter
 c_m=c{i}-mean(c{i});
 X=fft(c{i}.*winh,15000);
 bpm=abs(X);
 bpm(1:30)=0;            
 bpm(54:end)=0;
 [yvalue,k]=max(bpm);
 hr(i)=(k-1)*60*Fs/(length(X)); % Vector de ritmo cardiaco en BPM
 Per = round(15000/(k-1));
 %extracci�n de arm�nicos
 cSF=abs(fft(c_m(1:Per)))/yvalue;
 ARM(i,:)=cSF(3:6);
 end 
%%
%Time domain features for BVP
%Emotion
[amplitud,muestra,ancho,prominence]=findpeaks(pulsoaf,'MinPeakDistance',250,'MinPeakProminence',0.1,'MaxPeakWidth',250);
NN=diff(muestra)/500;
interbeat=60./NN;
%Baseline
[amplitud,muestra,ancho,prominence]=findpeaks(pulsolb,'MinPeakDistance',250,'MinPeakProminence',0.1,'MaxPeakWidth',250);
NNlb=diff(muestra)/500;
interbeatlb=60./NNlb;
%%%%
%THD Baseline Compute
thd2lb=mean(ARMlb(:,1));
thd3lb=mean(ARMlb(:,2));
thd4lb=mean(ARMlb(:,3));
thd5lb=mean(ARMlb(:,4));
%THD Emotion Compute
thd2=mean(ARM(:,1));
thd3=mean(ARM(:,2));
thd4=mean(ARM(:,3));
thd5=mean(ARM(:,4));

%%
%GSR time domain features

diff_gsrmm_lb=diff(gsrmmlb); %Derivate of the baseline GSR
avg_diff_gsrmm_lb=mean(diff_gsrmm_lb); %Average of the GSR's baseline derivate
amountneg_diff_gsrmm_lb=length(diff_gsrmm_lb(diff_gsrmm_lb<0));%Amount of negative samples
ratio_amountneg_lb=amountneg_diff_gsrmm_lb/length(diff_gsrmm_lb);%Proportion of negative values vs the whole lb values
%------------------------------------------------------------------------------------------
diff_gsrmm   =diff(gsrmm);   %Derivate of the GSR during the stimulus
avg_dif_gsrmm = mean(diff_gsrmm); %%Average of the GSR response during the stimulus
amountneg_diff_gsrmm=length(diff_gsrmm(diff_gsrmm<0));%Amount of negative samples
ratio_amountneg=amountneg_diff_gsrmm/length(diff_gsrmm);%Proportion of negative values vs the whole stimulus values
% 
% GSR frequency domain features
% Computing IMFs for Baseline data
full_emd_lb=emd(gsrmmlb);
imf1.lb=full_emd_lb(1,:);
imf2.lb=full_emd_lb(2,:);
imf3.lb=full_emd_lb(3,:);
imf4.lb=full_emd_lb(4,:);
res.emd.lb=sum(full_emd_lb(5:end,:));
% 
% Energy of gsr signal
nmvgsrlb=gsrmmlb-mean(gsrmmlb);
e.signal=(sum((nmvgsrlb).^2))/length(nmvgsrlb);
% Energy of IMFs
e.imf1=(sum((imf1.lb).^2))/length(imf1.lb);
e.imf2=(sum((imf2.lb).^2))/length(imf2.lb);
e.imf3=(sum((imf3.lb).^2))/length(imf3.lb);
e.imf4=(sum((imf4.lb).^2))/length(imf4.lb);
e.res=(sum((res.emd.lb).^2))/length(res.emd.lb);
e.whole.imf=[e.imf1,e.imf2,e.imf3,e.imf4];
% zero crossing rate for
crmlb1=ZCR(imf1.lb);
crmlb2=ZCR(imf2.lb);
crmlb3=ZCR(imf3.lb);
crmlb4=ZCR(imf4.lb);

% %
% Computing IMFs for Emotion
% 
full_emd=emd(gsrmm);
imf1=full_emd(1,:);
imf2=full_emd(2,:);
imf3=full_emd(3,:);
imf4=full_emd(4,:);
res.emd=sum(full_emd(5:end,:));
% 
% Energy of gsr signal
nmvgsr=gsrmm-mean(gsrmm);
e.signalem=(sum((nmvgsr).^2))/length(nmvgsr);
% Energy of IMFs
e.imf1e=(sum((imf1).^2))/length(imf1);
e.imf2e=(sum((imf2).^2))/length(imf2);
e.imf3e=(sum((imf3).^2))/length(imf3);
e.imf4e=(sum((imf4).^2))/length(imf4);
e.rese=(sum((res.emd).^2))/length(res.emd);
emotion.whole.imf=[e.imf1e,e.imf2e,e.imf3e,e.imf4e];
% zero crossing rate for
crm1=ZCR(imf1);
crm2=ZCR(imf2);
crm3=ZCR(imf3);
crm4=ZCR(imf4);
%%
%Computing LF and HF HRV
[Pxxh,Fh] = periodogram(pulsolb,hann(length(pulsolb)),length(pulsolb),Fs);
ptoth = bandpower(Pxxh,Fh,'psd');
pbandHFh = bandpower(Pxxh,Fh,[0.16 0.4],'psd');
pbandLFh=  bandpower(Pxxh,Fh,[0.04 0.15],'psd');
ratioh=pbandLFh/pbandHFh;
%Normalized values from https://www.frontiersin.org/articles/10.3389/fphys.2014.00177/full
LFnu=1/(1+(ratioh)^-1);
HFnu=1/(1+(ratioh));

[Pxxe,Fe] = periodogram(pulsoaf,hann(length(pulsoaf)),length(pulsoaf),Fs);
ptote = bandpower(Pxxe,Fe,'psd');
pbandHFe= bandpower(Pxxe,Fe,[0.16 0.4],'psd');
pbandLFe=  bandpower(Pxxe,Fe,[0.04 0.15],'psd');
ratioe=pbandLFe/pbandHFe;
%Normalized values from https://www.frontiersin.org/articles/10.3389/fphys.2014.00177/full
LFnue=1/(1+(ratioe)^-1);
HFnue=1/(1+(ratioe));


disp('Features extraction')

disp('Baseline data')
 
%Array info

infolb=[mean(gsrmmlb),std(gsrmmlb),max(gsrmmlb)-min(gsrmmlb),avg_diff_gsrmm_lb,...
    amountneg_diff_gsrmm_lb,ratio_amountneg_lb,e.imf1,e.imf2,e.imf3,e.imf4,...
    crmlb1,crmlb2,crmlb3,crmlb4,mean(hrlb),std(hrlb),max(hrlb)-min(hrlb),...
    mode(hrlb),std(diff(NNlb)),rms(diff(NNlb)),thd2lb,thd3lb,thd4lb,thd5lb,LFnu,HFnu,ratioh];
% 
infoem=[mean(gsrmm),std(gsrmm),max(gsrmm)-min(gsrmm),avg_dif_gsrmm,...
    amountneg_diff_gsrmm,ratio_amountneg,e.imf1e,e.imf2e,e.imf3e,e.imf4e,...
    crm1,crm2,crm3,crm4,mean(hr),std(hr),max(hr)-min(hr),...
    mode(hr),std(diff(NN)),rms(diff(NN)),thd2,thd3,thd4,thd5,LFnue,HFnue,ratioe];

%By hand, average power in time domain
% f2=norm(pulsolb,2)^2/numel(pulsolb)



% 
 fprintf('%0.8f,',infolb)
 disp(',')
 fprintf('%0.8f,',infoem)


 
%%     
%Las señales rojas son las correspondientes a la prueba
% figure(1)
% subplot(2,2,1)
% plot(t1,gsrmm,'r')
% 
% subplot(2,2,2)
% plot(tglb,gsrmmlb,'b')
%  
% subplot(2,2,3);
% plot(t,pulsoaf,'r')
% 
% subplot(2,2,4); 
% plot(tplb,pulsolb,'b')
% 
% figure(2)
% subplot(2,1,1);
% plot(gsrmm,'r')
% subplot(2,1,2);
% plot(diff(gsrmm),'b')
%


% Baseline
% 
% GSR
% fprintf('SCR Neutral  Mean %8f .\n',mean(gsrmmlb));
% fprintf('SCR Neutral  Standard deviation %8f .\n',std(gsrmmlb));
% fprintf('SCR Neutral  Dynamic range %8f .\n',max(gsrmmlb)-min(gsrmmlb));
% fprintf('SCR Average Derivate %8f .\n',avg_diff_gsrmm_lb);
% fprintf('SCR Amount of negative samples %8f .\n',amountneg_diff_gsrmm_lb);
% fprintf('SCR Proportion of negative values vs the whole stimulus %8f .\n',ratio_amountneg_lb);
% 
% BVP
% fprintf('Heart Rate Neutral Mean %8f .\n',mean(hrlb));
% fprintf('Heart Rate Neutral Standard deviation %8f .\n',std(hrlb));
% fprintf('Heart Rate Neutral Dynamic range %8f .\n',(max(hrlb)-min(hrlb)));
% fprintf('Heart Rate Neutral Mode (Frequency) %8f .\n',mode(hrlb));
% fprintf('Heart Rate Mean (Time) %8f .\n',mean(interbeatlb));
% fprintf('Heart Rate SDNN (Time) %8f .\n',std(diff(NNlb)));
% fprintf('Heart Rate RMSSD (Time) %8f .\n',rms(diff(NNlb)));
% 
% Emotion
% 
% disp('Emotion')
% fprintf('SCR Mean %8f .\n',mean(gsrmm));
% fprintf('SCR Standard deviation %8f .\n',std(gsrmm));
% fprintf('SCR Neutral  Dynamic range %8f .\n',(max(gsrmm)-min(gsrmm)));
% fprintf('SCR Average Derivate %8f .\n',avg_dif_gsrmm);
% fprintf('SCR Amount of negative samples %8f .\n',amountneg_diff_gsrmm);
% fprintf('SCR Proportion of negative values vs the whole stimulus values %8f .\n',ratio_amountneg);
% 
% 
% fprintf('Heart Rate Mean (FFT) %8f .\n',mean(hr));
% fprintf('Heart Rate Standard deviation %8f .\n',std(hr));
% fprintf('Heart Rate Dynamic range %8f .\n',(max(hr)-min(hr)));
% fprintf('Heart Rate Mode (Frequency) %8f .\n',mode(hr));
% fprintf('Heart Rate Mean (Time) %8f .\n',mean(interbeat));
% fprintf('Heart Rate SDNN (Time) %8f .\n',std(diff(NN)));
% fprintf('Heart Rate RMSSD (Time) %8f .\n',rms(diff(NN)));
% 
