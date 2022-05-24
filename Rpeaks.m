%% PAN TOMPKING -> GOLD STANDARD FOR ECG
% R peaks e segnale filtrato 
clear all
close all
clc 

%%
folder = 'C:\Users\feder\Desktop\Tesi\Data\Filtered ECG';
list = dir(folder);
list(1) = [];
list(1) = [];
N = length(list)
list(N) = [];
N = length(list)

fs_Ecg = 1024;
% wind = 10*1024;
wind = 30*1024;
addpath 'C:\Users\feder\Desktop\Tesi\Data\Filtered ECG'
%%
for i = 1:1
    FOLDER = fullfile(list(i).folder, list(i).name);
    file = dir(FOLDER);
    name = file.name;
    load(name)
%     dim = size(Ecg.Values,1)
    dim = size(ECG_filt,1)
    n_window = round(dim/wind);
%     qrs_amp = zeros(1,dim);
    qrs_i = zeros(1,dim);
    ecg_filt = zeros(dim,1);
    n_peaks = 0;
   
    
    peaksFORwindow = zeros(n_window,1);

     for j = 1:n_window-1
%         [qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(Ecg.Values((j-1)*wind+1:j*wind),fs_Ecg,0);
        [~,qrs_i_raw,delay]=pan_tompkin(ECG_filt((j-1)*wind+1:j*wind),fs_Ecg,0);
        n_peaks_new = size(qrs_i_raw,2);
%         n_ecg_new = size(ecg_h,1)
        if j == 1
            qrs_i(:,n_peaks+1:n_peaks+n_peaks_new) = qrs_i_raw;
        else
            qrs_i(:,n_peaks+1:n_peaks+n_peaks_new) = qrs_i_raw+(j-1)*wind ; 
        end 
%         qrs_amp(:,n_peaks+1:n_peaks+n_peaks_new) = qrs_amp_raw;
%         ecg_filt(n_ecg+1:n_ecg+n_ecg_new) = ecg_h;
        n_peaks = n_peaks + n_peaks_new;
%         n_ecg = n_ecg + n_ecg_new; 
        peaksFORwindow(j) = n_peaks;
     end 
    num = n_peaks;
%     num_ecg = n_ecg;
%     qrs_AMP = zeros(1,num);
    qrs_I = zeros(1,num);
%     ecg_FILT = zeros(num_ecg,1);
%     qrs_AMP = qrs_amp(1:num);
    qrs_I = qrs_i(1:num); 
%     ecg_FILT = ecg_filt(1:num_ecg);

    qrs_AMP = zeros((peaksFORwindow(end-1)),1);
    for i = 1:peaksFORwindow(end-1)
        qrs_AMP(i) = ECG_filt(qrs_I(i));
    end 

    % Una volta che ho qrs_I, posso calcoare direttamente HR al minuto ed
    % ogni 5 minuti per 24 ore
    [HR_min,HR_5min] = HR(qrs_I,fs_Ecg);
   % Save 
     name_PT = erase(name,"ECG_FILT-")
     save(['C:\Users\feder\Desktop\Tesi\Data\Pan-Tompkins\' 'PT-' name_PT],'qrs_I','qrs_AMP','peaksFORwindow','HR_5min','HR_min')
     %save(['PT-' name_PT],'qrs_I','qrs_AMP','peaksFORwindow')
end


%% IMPORTANTE PER FARE UN CONTROLLO VELOCE 
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