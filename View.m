% View data
clear all
close all
clc
%%
folderS = 'C:\Users\feder\Desktop\Tesi\Data\24h Signals'
fullFileNameS = fullfile(folderS, '2021-01-16 15.01.53.mat');
load(fullFileNameS)
%% Check dimension
Ore_SCG = (size(Acc.X_Acc,1))/(64*3600)
Ore_ECG = (size(Ecg.Values,1))/(1024*3600)  
%% Check Signal
folderACCF = 'C:\Users\feder\Desktop\Tesi\Data\Filtered Signals'
fullFileNameACCF = fullfile(folderACCF, 'FILT-2021-01-16 15.01.53.mat');
load(fullFileNameACCF)

figure()
a(1) = subplot(311), plot((Acc_filt.Var4./3600),Acc_filt.x_filt),title('Acc filtrata x')
a(2) = subplot(312), plot((Acc_filt.Var4./3600),Acc_filt.y_filt),title('Acc filtrata y')
a(3) = subplot(313), plot((Acc_filt.Var4./3600),Acc_filt.z_filt),title('Acc filtrata z')
linkaxes(a,'x')

figure()
a(1) = subplot(311), plot((Rot_filt.Var4./3600),Rot_filt.x_rotfilt),title('Rot filtrata x')
a(2) = subplot(312), plot((Rot_filt.Var4./3600),Rot_filt.y_rotfilt),title('Rot filtrata y')
a(3) = subplot(313), plot((Rot_filt.Var4./3600),Rot_filt.z_rotfilt),title('Rot filtrata z')
linkaxes(a,'x')



%% Check PT 
addpath 'C:\Users\feder\Desktop\Tesi\Data'\Pan-Tompkins\
addpath 'C:\Users\feder\Desktop\Tesi\Data\Filtered ECG'
addpath 'C:\Users\feder\Desktop\Tesi\Data'\'Picchi SCG - Acc z'\
load 'PT-2022-04-08 12.15.40.mat'
load 'ECG_FILT-2022-04-08 12.15.40.mat'
load 'SCG(Az)_picchi-2021-05-20 20.00.00.mat'
% figure()
% plot(ecg_FILT)
% hold on 
% scatter(qrs_I,qrs_AMP,'m')

wind = 30*1024;
for i = 1:10
    if i == 1
        figure()
        plot((i-1)*wind+1:i*wind,ECG_filt((i-1)*wind+1:i*wind))
        hold on
        scatter(qrs_I(peaksFORwindow(i):peaksFORwindow(i)),qrs_AMP(peaksFORwindow(i):peaksFORwindow(i)),'m')
    else 
        figure()
        plot((i-1)*wind+1:i*wind,ECG_filt((i-1)*wind+1:i*wind))
        hold on
        scatter(qrs_I(peaksFORwindow(i-1)+1:peaksFORwindow(i)),qrs_AMP(peaksFORwindow(i-1)+1:peaksFORwindow(i)),'m')
    end 
end 
%%
% Quanti battiti ci sono al minuto (60 sec) e quanto in ogni finestra di 10 s?
n_finestre = size(peaksFORwindow,1);
n_picchi_min = zeros(round(n_finestre/2),1);
n_picchi_30s = zeros(n_finestre,1);
for i = 1:n_finestre-1 
    if i == 1
        n_picchi_30s(i) = peaksFORwindow(i);
    else 
        n_picchi_30s(i) = peaksFORwindow(i)-peaksFORwindow(i-1);
    end
end 

mean_picchi_30s = mean(n_picchi_30s)
mean_picchi_min = 2*mean_picchi_30s    







%%
% load ("FILT-2021-06-12 09.38.41 MAIOLO.mat")
% n = length(Acc_filt.X_Acc)
% 
% figure()
% subplot(321),plot(Acc.Time_Acc./3600,Acc.X_Acc),xlabel('Hour [H]'),title('Acceleration x')
% subplot(322),plot(Acc.Time_Acc./3600,Acc_filt.x_filt),xlabel('Hour [H]'),title('Filtered Acceleration x')
% subplot(323),plot(Acc.Time_Acc./3600,Acc.Y_Acc),xlabel('Hour [H]'),title('Acceleration y')
% subplot(324),plot(Acc.Time_Acc./3600,Acc_filt.y_filt),xlabel('Hour [H]'),title('Filtered Acceleration y') 
% subplot(325),plot(Acc.Time_Acc./3600,Acc.Z_Acc),xlabel('Hour [H]'),title('Acceleration z')
% subplot(326),plot(Acc.Time_Acc./3600,Acc_filt.z_filt),xlabel('Hour [H]'),title('Filtered Acceleration z') 
% sgtitle('Acceleration VS Filtered Acceleration')
% 
% figure()
% subplot(321),plot(Rot.Time_Rot./3600,Rot.X_Rot),xlabel('Hour [H]'),title('Angular Velocity x')
% subplot(322),plot(Rot.Time_Rot./3600,Rot_filt.x_rotfilt),xlabel('Hour [H]'),title('Filtered Angular Velocity x')
% subplot(323),plot(Rot.Time_Rot./3600,Rot.Y_Rot),xlabel('Hour [H]'),title('Angular Velocity y')
% subplot(324),plot(Rot.Time_Rot./3600,Rot_filt.y_rotfilt),xlabel('Hour [H]'),title('Filtered Angular Velocity y') 
% subplot(325),plot(Rot.Time_Rot./3600,Rot.Z_Rot),xlabel('Hour [H]'),title('Angular Velocity z')
% subplot(326),plot(Rot.Time_Rot./3600,Rot_filt.z_rotfilt),xlabel('Hour [H]'),title('Filtered Angular Velocity z') 
% sgtitle('Angular Velocity VS Filtered Angular Velocity')

figure()
subplot(311),plot(Acc.Time_Acc./3600,Acc.X_Acc),xlabel('Hour [H]'),ylabel('[millig]')
hold on; plot(Acc.Time_Acc./3600,Acc_filt.x_filt),legend('Acceleration X','Filtered Acceleration X'),
title('Acceleration - X axis')
subplot(312),plot(Acc.Time_Acc./3600,Acc.Y_Acc),xlabel('Hour [H]'),ylabel('[millig]')
hold on; plot(Acc.Time_Acc./3600,Acc_filt.y_filt),legend('Acceleration Y','Filtered Acceleration Y')
title('Acceleration - Y axis')
subplot(313),plot(Acc.Time_Acc./3600,Acc.Z_Acc),xlabel('Hour [H]'),ylabel('[millig]')
hold on; plot(Acc.Time_Acc./3600,Acc_filt.z_filt),legend('Acceleration Z','Filtered Acceleration Z')
title('Acceleration - Z axis')



%% CONTROLLO SCG - DEVO SCEGLIERE LA FINESTRA TEMPORALE GIUSTA PER IL CALCOLO DI SCG 
addpath 'C:\Users\feder\Desktop\Tesi\Data'\Pan-Tompkins\
addpath 'C:\Users\feder\Desktop\Tesi\Data\Filtered ECG'
addpath 'C:\Users\feder\Desktop\Tesi\Data'\'Picchi SCG - Acc z'\
% addpath 'C:\Users\feder\Desktop\Tesi\Data'\'Cambia'\
%%
load 'PT-2022-02-21 20.10.53.mat'
load 'ECG_FILT-2022-02-21 20.10.53.mat'
load 'SCG(Az)_picchi-2022-02-21 20.10.53.mat'
battiti_ECG = peaksFORwindow(end-1)
battici_ACCZ = peaksFORwindow_SCG30(end-1)
%%
wind_ecg = 30*1024;
wind_ke = 30*64;

n_finestre = 600;
figure()
a = subplot(211) 
plot(((n_finestre-1)*wind_ecg+1:(n_finestre)*wind_ecg)./(1*1024),ECG_filt((n_finestre-1)*wind_ecg+1:(n_finestre)*wind_ecg)),xlabel('[s]'),title('Ecg PT')
hold on
scatter((qrs_I(peaksFORwindow(n_finestre-1)+1:peaksFORwindow(n_finestre)))./(1*1024),qrs_AMP(peaksFORwindow(n_finestre-1)+1:peaksFORwindow(n_finestre)),'m')
b = subplot(212)
plot(((n_finestre-1)*wind_ke+1:n_finestre*wind_ke)./(1*64),Acc_z((n_finestre-1)*wind_ke+1:n_finestre*wind_ke)),xlabel('[s]'),title('Acc z filtrata - template matching')
hold on
plot((POS_picchi_SCG(peaksFORwindow_SCG30(n_finestre-1)+1:peaksFORwindow_SCG30(n_finestre)))./(64*1),AMP_picchi_SCG(peaksFORwindow_SCG30(n_finestre-1)+1:peaksFORwindow_SCG30(n_finestre)),'*r')
sgtitle('Finestra di 30 sec di notte - template selezionato di 10 sec')
linkaxes([a b],'x')






