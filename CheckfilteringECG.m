%% Check segnale ECG filtrato e non
 clear all
 close all
 clc
%% 
folderECGF = 'C:\Users\feder\Desktop\Tesi\Data\Filtered ECG';
folderECG = 'C:\Users\feder\Desktop\Tesi\Data\ECG';
folderPT = 'C:\Users\feder\Desktop\Tesi\Data\PostProc PT'; 

listECGF = dir(folderECGF);
listECGF(1) = [];
listECGF(1) = [];
listECGF(end) = [];
listECG = dir(folderECG);
listECG(1) = [];
listECG(1) = [];
listECG(end) = [];
listPT = dir(folderPT);
listPT(1) = [];
listPT(1) = [];
N = length(listECG);

addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Filtered ECG'\
addpath 'C:\Users\feder\Desktop\Tesi'\Data\'ECG'\
addpath 'C:\Users\feder\Desktop\Tesi'\Data\'PostProc PT'\
%%
for m = 9:9
    FOLDERECG = fullfile(listECG(m).folder, listECG(m).name)
    file = dir(FOLDERECG);
    name = file.name;
    load(name)

    FOLDERECGF = fullfile(listECGF(m).folder, listECGF(m).name)
    file = dir(FOLDERECGF);
    name = file.name;
    load(name)

    FOLDERPT = fullfile(listPT(m).folder, listPT(m).name)
    file = dir(FOLDERPT);
    name = file.name;
    load(name)

     % Aggiunto dopo il controllo di RR 
    qrs_I = R_post_processing(:,1)';
    qrs_AMP = R_post_processing(:,2);

    fs_ECG = 1024;
    fs_SCG = 64;
    window_SCG = 0.2*fs_SCG;
    window_ECG = 0.2*fs_ECG;

    iniziopicchi = sum(HR_min(1:5));
    finepicchi = sum(HR_min(end-5:end)); 
    picchi_totali = length(qrs_I)-iniziopicchi-finepicchi; 
    Ecg_valori = Ecg.Values;
%%
%     for i = iniziopicchi:length(qrs_I)-finepicchi-1
for i = 10000:10100
        finestrabattito_ECGF = ECG_filt(qrs_I(i)-window_ECG:qrs_I(i+1)-window_ECG)';
        finestrabattito_ECG = Ecg_valori(qrs_I(i)-window_ECG:qrs_I(i+1)-window_ECG)';
        figure()
        set(gcf, 'WindowState', 'maximized');
        a = subplot(211); plot((qrs_I(i)-window_ECG:qrs_I(i+1)-window_ECG)./1024,finestrabattito_ECG),xlabel('[s]'),title('Ecg')
        b = subplot(212); plot((qrs_I(i)-window_ECG:qrs_I(i+1)-window_ECG)./1024,finestrabattito_ECGF),xlabel('[s]'),title('Filtered Ecg')
        sgtitle(i)
        pause
        close all
    end

end


