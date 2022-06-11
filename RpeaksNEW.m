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
% N = length(list)
N = length(list)-1

fs_Ecg = 1024;
% wind = 10*1024;
wind = 30*1024;
addpath 'C:\Users\feder\Desktop\Tesi\Data\Filtered ECG'
%%
for i = 1:N
    FOLDER = fullfile(list(i).folder, list(i).name);
    file = dir(FOLDER);
    name = file.name;
    load(name)
    dim = size(ECG_filt,1)
    n_window = round(dim/wind);
    qrs_i = zeros(1,dim);
    ecg_filt = zeros(dim,1);
    n_peaks = 0;
    qrs_I = [];
    qrs_AMP = [];
    
    peaksFORwindow = zeros(n_window,1);

     for j = 1:n_window-1
        [~,qrs_i_raw,delay]=pan_tompkin(ECG_filt((j-1)*wind+1:j*wind),fs_Ecg,0);
        n_peaks_new = size(qrs_i_raw,2);
%         plot(ECG_filt((j-1)*wind+1:j*wind)); hold on; plot(qrs_i_raw,0,'*r')
%         n_ecg_new = size(ecg_h,1)
        if j == 1
            qrs_i(:,n_peaks+1:n_peaks+n_peaks_new) = qrs_i_raw;
        else
            qrs_i(:,n_peaks+1:n_peaks+n_peaks_new) = qrs_i_raw+(j-1)*wind ; 
        end 
     
       for p = n_peaks+1:n_peaks_new+n_peaks
            search = ECG_filt(qrs_i(p)-0.02*fs_Ecg:qrs_i(p)+0.02*fs_Ecg,1);
            new_max = max(search);
            pos_max = find(search == new_max);
            pos = qrs_i(p)-0.02*fs_Ecg+pos_max-0.88;

%             pos_max1 = find(ECG_filt == new_max); % VA BENE MA MEGLIO CERCARLO SOLO NELLA FINESTRA, potrei avere più punti nel segnale con lo stesso valore 277106
%             diff = pos-pos_max1;
%             figure() 
%             plot(qrs_i(p)-0.02*fs_Ecg:qrs_i(p)+0.02*fs_Ecg,search); hold on; plot(qrs_i(p),ECG_filt(qrs_i(p)),'*r'); hold on; plot(pos_max1,new_max,'*m'); hold on; plot(pos,new_max,'.b')
%             pause
%             close all


            qrs_AMP = [qrs_AMP; new_max];
            qrs_I = [qrs_I pos];
       end
        n_peaks = n_peaks + n_peaks_new;
        peaksFORwindow(j) = n_peaks;
     end 
    num = n_peaks;

%     qrs_I = zeros(1,num);
%     qrs_I = qrs_i(1:num); 
% 
%     qrs_AMP = zeros((peaksFORwindow(end-1)),1);
%     for i = 1:peaksFORwindow(end-1)
%         qrs_AMP(i) = ECG_filt(qrs_I(i));
%     end 


    [HR_min,HR_5min] = HR(qrs_I,fs_Ecg);
   % Save 
     name_PT = erase(name,"ECG_FILT-")
     save(['C:\Users\feder\Desktop\Tesi\Data\Pan-Tompkins\' 'PT-' name_PT],'qrs_I','qrs_AMP','peaksFORwindow','HR_5min','HR_min')
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
    pause
end 
