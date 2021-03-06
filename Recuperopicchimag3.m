%% Federica Mozzini - 946400
% Recupero i picchi che PT ha perso durante l'analisi
 clear all
 close all
 clc
%%
% folderSCG = 'C:\Users\feder\Desktop\Tesi\Data\Picchi SCG - Acc z';
% folderPT = 'C:\Users\feder\Desktop\Tesi\Data\PostProc PT'; 
folderECG = 'C:\Users\feder\Desktop\Tesi\Data\Filtered ECG';
folderSCG = 'C:\Users\feder\Desktop\Tesi\Data\Picchi SCG - Acc z - 10 sec';
folderPT = 'C:\Users\feder\Desktop\Tesi\Data\PostProc PT_10SEC'; 

listPT = dir(folderPT);
listPT(1) = [];
listPT(1) = [];
listECG = dir(folderECG);
listECG(1) = [];
listECG(1) = [];
listECG(end) = [];
list = dir(folderSCG);
list(1) = [];
list(1) = [];
N = length(list)-1;
% list(N-1) = [];
% list(N-1) = [];
% N = length(list)

 addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Filtered ECG'\
%  addpath 'C:\Users\feder\Desktop\Tesi'\Data\'PostProc PT'\
%  addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Picchi SCG - Acc z'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Codes\
 addpath 'C:\Users\feder\Desktop\Tesi'\Data\'PostProc PT_10SEC'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Picchi SCG - Acc z - 10 sec'\
 %%
 for m = 1:1
    FOLDERSCG = fullfile(list(m).folder, list(m).name)
    file = dir(FOLDERSCG);
    name = file.name;
    load(name)

    FOLDERPT = fullfile(listPT(m).folder, listPT(m).name)
    file = dir(FOLDERPT);
    name = file.name;
    load(name)

    FOLDERECG = fullfile(listECG(m).folder, listECG(m).name)
    file = dir(FOLDERECG);
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
    n_picchi = zeros(length(qrs_I)-finepicchi-1,1);
    qrs_new = [];
    error = 0;
%%
    for i = iniziopicchi:length(qrs_I)-finepicchi-1 %scorro tutti i picchi ECG   
        qrs1 = (qrs_I(i)/1024)*64;
        qrs2 = (qrs_I(i+1)/1024)*64;
        finestra_dopo = qrs2-window_SCG; %cambiato
        finestra_prima = qrs1-window_SCG;
        [row]=find(POS_picchi_SCG<finestra_dopo & POS_picchi_SCG>finestra_prima);
        n_picchi(i,1) = size(row,1);
        if n_picchi(i) >= 3
            finestrabattito_ECGF = ECG_filt(qrs_I(i)-window_ECG:qrs_I(i+1)-window_ECG)';
            finestrabattito_SCG = Acc_z(qrs1-window_SCG:qrs2-window_SCG)';
            try 
            [~,qrs_i_raw,delay]=pan_tompkin(finestrabattito_ECGF,fs_ECG,0);
            n_peaks_new = size(qrs_i_raw,2);
            inizio = qrs_I(i)-qrs_i_raw(1);
            qrs_i_new = inizio+qrs_i_raw;
            qrs_amp_new = ECG_filt(qrs_i_new);
            qrs_new1 = [qrs_i_new' qrs_amp_new];
            qrs_new = [qrs_new; qrs_new1];
%             figure()
%             set(gcf, 'WindowState', 'maximized');
%             a = subplot(211); plot((R(i,1)-window_ECG:R(i+1,1)-window_ECG)./1024,finestrabattito_ECGF),xlabel('[s]'),title('Ecg'); hold on;
%             plot(R(i,1)/1024,R(i,2),'*r'); hold on;
%             for p = 1: length(qrs_i_new)
%                 plot(qrs_new1(p,1)/1024,qrs_new1(p,2),'*g')
%             end 
%             b = subplot(212); plot((qrs1-window_SCG:qrs2-window_SCG)./64,finestrabattito_SCG),xlabel('[s]'),title('SCG'); hold on;
%             for r = 1:length(row)
%                 plot((POS_picchi_SCG(row(r)))./64,AMP_picchi_SCG(row(r)),'mo')
%             end
%             sgtitle(i)
%             pause 
%             close all
            catch ME 
                error = error+1;
            end
        end 
    end 
    qrs = [qrs_I' qrs_AMP];
    QRS = [qrs; qrs_new];
    QRS = sortrows(QRS);
    [HR_min,~,~] = HR(QRS(:,1),fs_ECG);
    RR = diff(QRS(:,1)); %campioni
    RR_s = RR./fs_ECG; %sec
    HR = 60./RR_s; %espressi in sec ma sono campioni ( ho 1 valore per ogni campione)
    HR_30heartbeat = zeros(round(length(HR)/30),1);
%         HR_30sec = zeros(length(HR),1);
    m = 0;
    for i = 1:30:length(HR)-30
        m = m+1;
        HR_30heartbeat(m) = mean(HR(i:i+29));
    end 


%     % Save
    name = erase(name,"ECG_FILT-")
% % %     save(['C:\Users\feder\Desktop\Tesi\Data\PostProc PT 1\' 'PostProc PT 1-' name], 'QRS','HR_min','HR_5min')
    save(['C:\Users\feder\Desktop\Tesi\Data\PostProc PT 1_10SEC\' 'PostProc PT 1-' name], 'QRS','HR_30heartbeat','HR_min')
    
 end 