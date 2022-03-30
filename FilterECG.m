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
for i = 1:N
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
wind_ecg = 10*60*1024;
for i = 30:36
figure()
plot(((i-1)*wind_ecg+1:(i)*wind_ecg)./(1*1024),ECG_filt((i-1)*wind_ecg+1:(i)*wind_ecg)),xlabel('[s]'),title('Ecg PT')
end