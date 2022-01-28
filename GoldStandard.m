%% Compute GOLD Standard
clear all
close all
clc

%% 
load 'RotationalKE-2021-01-16 15.01.53.mat'
load 'PT-2021-01-16 15.01.53.mat'

% 1O Hour -> 3600 secondi in 1 ora/10 s/w = 360 wind all'ora 
% alla decima ora ho 10*360
windPERhour = 360;
wind_ecg = 10*1024;
wind_sig = 10*64;
% Quante finestre ho in 10 ore? 

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
plot(cineticar(n_fin_iniziale*wind_sig:n_fin_finale*wind_sig)),title('Rotational Kinetic Energy 2:00 per 10 sec')
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