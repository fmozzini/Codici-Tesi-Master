%% Compute GOLD Standard
clear all
close all
clc

%% Qui devo farlo dentro ad un ciclo che va in entrambe le cartelle 
% Load Rotational KE
folderKE = 'C:\Users\feder\Desktop\Tesi\Data\RotationalKE'
fullFileNameKE = fullfile(folderKE, 'RotationalKE-2021-01-16 15.01.53.mat');
load(fullFileNameKE)

% Load Pan-Tompkin
folderPT = 'C:\Users\feder\Desktop\Tesi\Data\Pan-Tompkins'
fullFileNamePT = fullfile(folderPT, 'PT-2021-01-16 15.01.53.mat');
load(fullFileNamePT)

% Load Filtered Acc per confronto
folderACCF = 'C:\Users\feder\Desktop\Tesi\Data\Filtered Signals'
fullFileNameACCF = fullfile(folderACCF, 'FILT-2021-01-16 15.01.53.mat');
load(fullFileNameACCF)
%% A partire dalla KE devo cercare di trovare i picchi - RKE
dim = size(cineticar,1); %durata in campioni
fs_KE = 64;
wind_30sec = fs_KE*30;
n_window = round(dim/wind_30sec);
% template = [];
HR_picchi = zeros(dim,1);
Pos_picchi = zeros(dim,1);
Amp_picchi = zeros(dim,1);
n_peaks = 0;
peaksFORwindow_RKE30 = zeros(n_window,1); 
% Identifico il template di 10 sec durante la notte o molto più corto solo
% del segnale che voglio
figure()
wind_10sec = 64*10;
plot(((4000-1)*wind_10sec+1:4000*wind_10sec)./(1*64),cineticar((4000-1)*wind_10sec+1:4000*wind_10sec)),title('Template 10 sec - RKE'),xlabel('[s]')
brush('on')
figure()
plot(brushedData)
templateBR = brushedData(:,2);
templateCIN = cineticar((4000-1)*wind_10sec+1:4000*wind_10sec);


subplot(211)
plot(templateBR),title('Template selezionato di poco più di 1 sec')
subplot(212)
plot(templateCIN),title('10 sec')

%%
for i = 1:n_window-1 
% for i = 1:2 
    [HR, pos_picchi, amp_picchi,Template] = SCG_template_matching_corr(cineticar((i-1)*wind_30sec+1:wind_30sec*i),fs_KE,'x',templateBR)
    n_peaks_new = size(pos_picchi,1);
    if i == 1
            Pos_picchi(n_peaks+1:n_peaks+n_peaks_new,:) = pos_picchi;
        else
            Pos_picchi(n_peaks+1:n_peaks+n_peaks_new,:) = pos_picchi+(i-1)*wind_30sec ; 
        end 
    Amp_picchi(n_peaks+1:n_peaks+n_peaks_new,:) = amp_picchi;
%     HR_picchi(n_peaks+1:n_peaks+n_peaks_new,:) = HR;
    n_peaks = n_peaks + n_peaks_new;
    peaksFORwindow_RKE30(i) = n_peaks;
end 
num = n_peaks;
AMP_picchi_RKE = zeros(num,1);
POS_picchi_RKE = zeros(num,1);
%HR_picchi_RKE = zeros(num,1);
AMP_picchi_RKE = Amp_picchi(1:num);
POS_picchi_RKE = Pos_picchi(1:num); 
%HR_picchi_RKE = HR_picchi(1:num);

%% PROVO SU ACC Z
dim1 = size(Acc_filt.z_filt,1); %durata in campioni
fs_KE1 = 64;
wind_30sec1 = fs_KE1*30;
n_window1 = round(dim1/wind_30sec1);
% template = [];
HR_picchi1 = zeros(dim1,1);
Pos_picchi1 = zeros(dim1,1);
Amp_picchi1 = zeros(dim1,1);
n_peaks1 = 0;
peaksFORwindow_RKE301 = zeros(n_window1,1); 
% Identifico il template di 10 sec durante la notte 
wind_10sec1 = 64*10;
figure()
plot(((4000-1)*wind_10sec1+1:4000*wind_10sec1)./(60*64),Acc_filt.z_filt((4000-1)*wind_10sec1+1:4000*wind_10sec1)),title('Template 10 sec - SCG:Acc filtrata z')
template1 = Acc_filt.z_filt((4000-1)*wind_10sec1+1:4000*wind_10sec1);


for i = 1:n_window1-1 
% for i = 1:2 
    [HR1, pos_picchi1, amp_picchi1,Template1] = SCG_template_matching_corr(Acc_filt.z_filt((i-1)*wind_30sec1+1:wind_30sec1*i),fs_KE1,'x',template1)
    n_peaks_new1 = size(pos_picchi1,1);
    if i == 1
            Pos_picchi1(n_peaks1+1:n_peaks1+n_peaks_new1,:) = pos_picchi1;
        else
            Pos_picchi1(n_peaks1+1:n_peaks1+n_peaks_new1,:) = pos_picchi1+(i-1)*wind_30sec1 ; 
        end 
    Amp_picchi1(n_peaks1+1:n_peaks1+n_peaks_new1,:) = amp_picchi1;
%     HR_picchi(n_peaks+1:n_peaks+n_peaks_new,:) = HR;
    n_peaks1 = n_peaks1 + n_peaks_new1;
    peaksFORwindow_RKE301(i) = n_peaks1;
end 
num1 = n_peaks1;
AMP_picchi_RKE1 = zeros(num1,1);
POS_picchi_RKE1 = zeros(num1,1);
%HR_picchi_RKE = zeros(num,1);
AMP_picchi_RKE1 = Amp_picchi1(1:num1);
POS_picchi_RKE1 = Pos_picchi1(1:num1); 
%HR_picchi_RKE = HR_picchi(1:num);
%% 
% PER CONFRONTARLI - POTREI FARE UN SUBPLOT
for i = 1:10
    if i == 1
        figure()
        plot((i-1)*wind+1:i*wind,ecg_FILT((i-1)*wind+1:i*wind))
        hold on
        scatter(qrs_I(peaksFORwindow(i):peaksFORwindow(i)),qrs_AMP(peaksFORwindow(i):peaksFORwindow(i)),'m')
    else 
        figure()
        plot((i-1)*wind+1:i*wind,ecg_FILT((i-1)*wind+1:i*wind))
        hold on
        scatter(qrs_I(peaksFORwindow(i-1)+1:peaksFORwindow(i)),qrs_AMP(peaksFORwindow(i-1)+1:peaksFORwindow(i)),'m')
    end 
end


%% Inizio Holter 15:00
% Notte (2-3 per essere sicuri) - 11° 12° ora 

% 11°esima ora, uso finestre di 30 sec
% 11*3600 = 39600/60 = 1320
% 39600/10 = 3960
wind_ecg = 30*1024;
wind_ke = 30*64;

figure()
subplot(311)
plot(((1320-1)*wind_ecg+1:(1320)*wind_ecg)./(1*1024),ecg_FILT((1320-1)*wind_ecg+1:(1320)*wind_ecg)),xlabel('[s]'),title('Ecg filtrato PT')
hold on
scatter((qrs_I(peaksFORwindow(1320-1)+1:peaksFORwindow(1320)))./(1*1024),qrs_AMP(peaksFORwindow(1320-1)+1:peaksFORwindow(1320)),'m')
subplot(312)
plot(((1320-1)*wind_ke+1:1320*wind_ke)./(1*64),cineticar((1320-1)*wind_ke+1:1320*wind_ke)),xlabel('[s]'),title('Rotational KE - template matching con 0,6*fs')
hold on
plot((POS_picchi_RKE(peaksFORwindow_RKE30(1320-1)+1:peaksFORwindow_RKE30(1320)))./(64*1),AMP_picchi_RKE(peaksFORwindow_RKE30(1320-1)+1:peaksFORwindow_RKE30(1320)),'*r')
subplot(313)
plot(((1320-1)*wind_ke+1:1320*wind_ke)./(1*64),Acc_filt.z_filt((1320-1)*wind_ke+1:1320*wind_ke)),xlabel('[s]'),title('Acc z filtrata - template matching')
hold on
plot((POS_picchi_RKE1(peaksFORwindow_RKE301(1320-1)+1:peaksFORwindow_RKE301(1320)))./(64*1),AMP_picchi_RKE1(peaksFORwindow_RKE301(1320-1)+1:peaksFORwindow_RKE301(1320)),'*r')
if isempty(template)
    sgtitle('Finestra di 30 sec di notte - No template')
else 
    sgtitle('Finestra di 30 sec di notte - template selezionato di poco più di 1 s')
end 

%%
% Hour = 3600*64;
% Ecg_Hour = 3600*1024;
% figure()
% subplot(211)
% plot(cineticar(10*Hour:10*Hour+630)),title('Rotational Kinetic Energy 2:00 per 10 sec')
% subplot(212)
% plot(ecg_FILT(10*Ecg_Hour:10*Ecg_Hour+10240)),title('Ecg 2:00 per 10 sec')

%% Faccio solo un ciclo che cicla le finestre
wind_ecg = 10*1024;
wind_sig = 10*64;
n_fin_iniziale = 3600;
n_fin_finale = 3601;

figure()
subplot(211)
plot(n_fin_iniziale*wind_sig:n_fin_finale*wind_sig,cineticar(n_fin_iniziale*wind_sig:n_fin_finale*wind_sig)),title('Rotational Kinetic Energy 2:00 per 10 sec')
subplot(212)
plot(n_fin_iniziale*wind_ecg:n_fin_finale*wind_ecg,ecg_FILT(n_fin_iniziale*wind_ecg:n_fin_finale*wind_ecg)),title('ECG and R peaks per 10 sec')
hold on 
scatter(qrs_I(peaksFORwindow(n_fin_iniziale):peaksFORwindow(n_fin_finale)),qrs_AMP(peaksFORwindow(n_fin_iniziale):peaksFORwindow(n_fin_finale)),'m')
%% FUNZIONA MA SI PU' FARE MEGLIO cd 
% wind_ecg = 10*1024;
% wind_sig = 10*64;
% finFORhour = 360;
% 
% figure()
% subplot(211)
% plot(cineticar(10*finFORhour*wind_sig:10*finFORhour*wind_sig+wind_sig)),title('Rotational Kinetic Energy 2:00 per 10 sec')
% subplot(212)
% plot(10*finFORhour*wind_ecg:10*finFORhour*wind_ecg+wind_ecg,ecg_FILT(10*finFORhour*wind_ecg:10*finFORhour*wind_ecg+wind_ecg)),title('ECG and R peaks per 10 sec')
% hold on 
% scatter(qrs_I(peaksFORwindow(10*finFORhour):peaksFORwindow(10*finFORhour+1)),qrs_AMP(peaksFORwindow(10*finFORhour):peaksFORwindow(10*finFORhour+1)),'m')