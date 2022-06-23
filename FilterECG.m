%% Filtered ECG (remove wander baseline)
clear all
close all 
clc

%%
folder = 'C:\Users\feder\Desktop\Tesi\Data\ECG';
list = dir(folder);
list(1) = [];
list(1) = [];
N = length(list);
list(N) = [];
%list(N-1) = [];
N = length(list)
addpath 'C:\Users\feder\Desktop\Tesi\Data\ECG'
%%
for i = N:N
    FOLDER = fullfile(list(i).folder, list(i).name);
    file = dir(FOLDER);
    name = file.name;
    load(name)

    Fs = 1024;
    ECG_filt = filtECG(Ecg.Values,Fs);

    name_ECG = erase(name,"ECG-")
    save(['C:\Users\feder\Desktop\Tesi\Data\Filtered ECG\' 'ECG_FILT-' name_ECG ],'ECG_filt')
end 

%% 
addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Filtered ECG'\
addpath 'C:\Users\feder\Desktop\Tesi'\Data\'ECG'\
load 'ECG_FILT-2021-01-16 15.01.53.mat'
load 'ECG-2021-01-16 15.01.53.mat'

figure()
subplot(121), plot(Ecg.Values)
subplot(122), plot(ECG_filt)
%%
wind_ecg = 10*1024;
for i = 405:410
figure()
a = subplot(211),plot(((i-1)*wind_ecg+1:(i)*wind_ecg)./(1*1024),Ecg.Values((i-1)*wind_ecg+1:(i)*wind_ecg)),xlabel('[s]'),ylabel('[mV]'); hold on;
b = subplot(212),plot(((i-1)*wind_ecg+1:(i)*wind_ecg)./(1*1024),ECG_filt((i-1)*wind_ecg+1:(i)*wind_ecg)),xlabel('[s]'),ylabel('[mV]'),title('Filtered Ecg')
linkaxes([a b],'x')
pause
end
%%
for i = 405:410
figure()
plot(((i-1)*wind_ecg+1:(i)*wind_ecg)./(1*1024),Ecg.Values((i-1)*wind_ecg+1:(i)*wind_ecg)),xlabel('[s]'); hold on;
plot(((i-1)*wind_ecg+1:(i)*wind_ecg)./(1*1024),ECG_filt((i-1)*wind_ecg+1:(i)*wind_ecg)),xlabel('[s]'),ylabel('[mV]'); legend('Ecg','Filtered Ecg'),
title('Ecg and Filtered Ecg')
pause
end