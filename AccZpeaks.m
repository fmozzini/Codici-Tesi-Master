%% SCG Peaks - Z axis filtered acceleration

clear all
close all
clc

%%
folderSCG = 'C:\Users\feder\Desktop\Tesi\Data\Filtered Signals';

list = dir(folderSCG);
list(1) = [];
list(1) = [];
N = length(list);
list(N-1) = [];
list(N-1) = [];
N = length(list)

[~,txtdata] = xlsread('C:\Users\feder\Desktop\Tesi\Info Pazienti.xlsx','H:I');
txtdata(1,:) = [];
Inizio_Holter = txtdata(:,1);
Periodo_Sonno = txtdata(:,2);
addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Filtered Signals'\
% addpath 'C:\Users\feder\Desktop\Tesi'\Data\'PostProc SCG'\
addpath 'C:\Users\feder\Desktop\Tesi'\Codes\
%% 
for m = 1:1
    FOLDERSCG = fullfile(list(m).folder, list(m).name)
    file = dir(FOLDERSCG);
    name = file.name;
    load(name)


    dim = size(Acc_filt.z_filt,1); %durata in campioni
    fs_SCG = 64;
    wind_30sec = fs_SCG*30;
    n_window = round(dim/wind_30sec);
    % template = [];
    HR_picchi = zeros(dim,1);
    Pos_picchi = zeros(dim,1);
    Amp_picchi = zeros(dim,1);
    n_peaks = 0;
    peaksFORwindow_SCG30 = zeros(n_window,1); 

    % Template
    % Cerco la notte 
    Hour_add = 2*3600;
    Inizio = datevec(Inizio_Holter{m}); 
    Sonno = datevec(Periodo_Sonno{m}); 
    durata_sec = etime(Sonno,Inizio)
    durata_h = durata_sec/3600
    Night = durata_h + 2;
    Night_sec = round(Night*3600*64);
%     templateNight = Acc_filt.z_filt(Night_sec:Night_sec+30*fs_SCG);
%     plot(Night_sec:Night_sec+30*fs_SCG,templateNight),title('Template di 30 sec preso durante la notte')
    for n = 1:2:200
        templateNight = Acc_filt.z_filt(Night_sec+30*n*fs_SCG:Night_sec+30*(n+1)*fs_SCG);
        plot(Night_sec+30*n*fs_SCG:Night_sec+30*(n+1)*fs_SCG,templateNight),title('Template di 30 sec preso durante la notte')
        pause
        answer = questdlg('Do you want to use this template?', ...
	        'SCG template matching template ', ...
	        'Yes','No','Yes');
        switch answer
            case 'No'
                continue
            case 'Yes'
                TemplateNight = templateNight;
                break 
        end 
    end 
    [~,~,~,Template] = SCG_template_matching_corr(TemplateNight,fs_SCG,'x',[]);
    plot(Template)
    pause
    for i = 1:n_window-1 
        [~, pos_picchi, amp_picchi,~] = SCG_template_matching_corr(Acc_filt.z_filt((i-1)*wind_30sec+1:wind_30sec*i),fs_SCG,'x',Template);
        pause
        n_peaks_new = size(pos_picchi,1);
        if i == 1
            Pos_picchi(n_peaks+1:n_peaks+n_peaks_new,:) = pos_picchi;
        else
            Pos_picchi(n_peaks+1:n_peaks+n_peaks_new,:) = pos_picchi+(i-1)*wind_30sec ; 
        end 
        Amp_picchi(n_peaks+1:n_peaks+n_peaks_new,:) = amp_picchi;
        n_peaks = n_peaks + n_peaks_new;
        peaksFORwindow_SCG30(i) = n_peaks;  
    end 
    num = n_peaks;
    AMP_picchi_SCG = zeros(num,1);
    POS_picchi_SCG = zeros(num,1);
    AMP_picchi_SCG = Amp_picchi(1:num);
    POS_picchi_SCG = Pos_picchi(1:num); 
    
    HR_SCG = 60./(diff(POS_picchi_SCG))*fs_SCG;

    % Save
%     name_SCG = erase(name,"FILT-")
%     Acc_z = Acc_filt.z_filt;
%     save(['C:\Users\feder\Desktop\Tesi\Data\Picchi SCG - Acc z\' 'SCG(Az)_picchi-' name_SCG],'Acc_z','AMP_picchi_SCG','POS_picchi_SCG','peaksFORwindow_SCG30','HR_SCG')
%     save(['C:\Users\feder\Desktop\Tesi\Data\Picchi SCG - Acc z\Template\' 'template-' name_SCG],'Template','templateNight')
end 

%% Numero di picchi SCG rilevati
battici_ACCZ = peaksFORwindow_SCG30(end-1)

%% Controllo se ?? tutto ok
folderPT = 'C:\Users\feder\Desktop\Tesi\Data\Pan-Tompkins'
fullFileNamePT = fullfile(folderPT, 'PT-2021-01-16 15.01.53.mat');
load(fullFileNamePT)
folderECG = 'C:\Users\feder\Desktop\Tesi\Data\Filtered ECG'
fullFileNameECG = fullfile(folderECG, 'ECG_FILT-2021-01-16 15.01.53.mat');
load(fullFileNameECG)
folderSCG = 'C:\Users\feder\Desktop\Tesi\Data\Picchi SCG - Acc z'
fullFileNameSCG = fullfile(folderSCG, 'SCG(Az)_picchi-2021-01-16 15.01.53.mat');
load(fullFileNameSCG)

wind_ecg = 30*1024;
wind_ke = 30*64;
%%
% Battiti ECG 
battiti_ECG = peaksFORwindow(end-1)
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

%% Check durante la notte
% n_finestre = 1000;
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
%% Check durante il giorno
% n_finestre = 600; % 5 ora
% n_finestre = 300; % 2 e 30 ora
n_finestre = 2040
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
%%
% load Pan-Tompkins\'PT-2021-06-24 16.30.32 MORRONE.mat'
% load 'SCG(Az)_picchi-2022-01-04 20.00.01 Riccardo Monti.mat'
% battiti_ECG = peaksFORwindow(end-1)
% battici_ACCZ = peaksFORwindow_SCG30(end-1)

