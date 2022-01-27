%% PAN TOMPKING -> GOLD STANDARD FOR ECG
% R peaks e segnale filtrato 
clear all
close all
clc 

%%
folder = 'C:\Users\feder\Desktop\Tesi\Data\ECG';
list = dir(folder);
list(1) = [];
list(1) = [];
list(14) = [];
N = length(list)

fs_Ecg = 1024;
wind = 10*1024;

for i = 1:1
    FOLDER = fullfile(list(1).folder, list(1).name);
    file = dir(FOLDER);
    name = file.name;
    load(name)
    dim = size(Ecg.Values,1)
    n_window = round(dim/wind);
    qrs_amp = zeros(1,dim);
    qrs_i = zeros(1,dim);
    ecg_filt = zeros(dim,1);
    n_peaks = 0;
    n_ecg = 0;
    
    peaksFORwindow = zeros(n_window,1);

    %for j = 1:n_window-1
       for j = 10:11
        [qrs_amp_raw,qrs_i_raw,delay,ecg_h]=pan_tompkin1(Ecg.Values((j-1)*wind+1:j*wind),fs_Ecg,0);
        n_peaks_new = size(qrs_i_raw,2);
        n_ecg_new = size(ecg_h,1)
        if j == 1
            qrs_i(:,n_peaks+1:n_peaks+n_peaks_new) = qrs_i_raw;
        else
            qrs_i(:,n_peaks+1:n_peaks+n_peaks_new) = qrs_i_raw+(j-1)*wind ; 
        end 
        qrs_amp(:,n_peaks+1:n_peaks+n_peaks_new) = qrs_amp_raw;
        ecg_filt(n_ecg+1:n_ecg+n_ecg_new) = ecg_h;
        n_peaks = n_peaks + n_peaks_new;
        n_ecg = n_ecg + n_ecg_new; 
        peaksFORwindow(j) = n_peaks;
     end 
    num = n_peaks;
    num_ecg = n_ecg;
    qrs_AMP = zeros(1,num);
    qrs_I = zeros(1,num);
    ecg_FILT = zeros(num_ecg,1);
    qrs_AMP = qrs_amp(1:num);
    qrs_I = qrs_i(1:num); 
    ecg_FILT = ecg_filt(1:num_ecg);

   % Save 
%     name_PT = erase(name,"ECG-")
%     save(['PT-' name_PT],'qrs_I','qrs_AMP',"ecg_FILT",'peaksFORwindow')
end

%%

figure()
plot(ecg_FILT)
hold on 
scatter(qrs_I,qrs_AMP,'m')
%%
peaksFORwindow(10) = 118 % numero di picchi
peaksFORwindow(20) = 265 % numero di picchi

figure()
plot(ecg_FILT(10*wind:20*wind))
hold on 
scatter(qrs_I(peaksFORwindow(10):peaksFORwindow(20)),qrs_AMP(peaksFORwindow(10):peaksFORwindow(20)),'m')