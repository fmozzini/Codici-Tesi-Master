%% Cut signals - Elimino i primi 5 minuti dei segnali, sono troppo corrotti

%% Prendo le finestre! Non considero i primi 5 minuti 
 clear all
 close all
 clc

folderSCG = 'C:\Users\feder\Desktop\Tesi\Data\Picchi SCG - Acc z';
folderPT = 'C:\Users\feder\Desktop\Tesi\Data\Pan-Tompkins';
folderECG = 'C:\Users\feder\Desktop\Tesi\Data\Filtered ECG';
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
N = length(list);
list(N-1) = [];
list(N-1) = [];
N = length(list)

 addpath 'C:\Users\feder\Desktop\Tesi'\Data\''\'Filtered ECG'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Data\''\Pan-Tompkins\
 addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Picchi SCG - Acc z'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Codes\


%% Parto da ECG, per ogni battito di ECG considero una finestra a sx e a dx del picco R.. 
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

    fs_ECG = 1024;
    fs_SCG = 64;
    window_SCG = 0.2*fs_SCG;
    window_ECG = 0.2*fs_ECG;

    %% non prendo dal primo picco - taglio i primi 5 minuti 
    cinqueminuti = 150;
    iniziopicchi = peaksFORwindow(cinqueminuti);
    
    FINESTREBATTITO_ECG = zeros(peaksFORwindow(end)-iniziopicchi,2000);
    FINESTREBATTITO_SCG = zeros(peaksFORwindow(end)-iniziopicchi,100);
 %%
 for i = iniziopicchi:peaksFORwindow(end-1)-1 %scorro tutti i picchi ECG

    qrs1 = (qrs_I(i)/1024)*64;
    qrs2 = (qrs_I(i+1)/1024)*64;
    finestra_dopo = qrs2;
    finestra_prima = qrs1-window_SCG;
    [row]=find(POS_picchi_SCG<finestra_dopo & POS_picchi_SCG>finestra_prima);
    n_picchi = size(row,1);
    
    picchi = POS_picchi_SCG(row);
    finestrabattito_ECG = ECG_filt(qrs_I(i)-window_ECG:qrs_I(i+1))';
    finestrabattito_SCG = Acc_z(qrs1-window_SCG:qrs2)';

    FINESTREBATTITO_ECG(i-iniziopicchi+1,1:length(finestrabattito_ECG)) = finestrabattito_ECG;
    FINESTREBATTITO_SCG(i-iniziopicchi+1,1:length(finestrabattito_SCG)) = finestrabattito_SCG;
%  
%     Opzione con cell array, non devo fare l'inizializzazione 
%     FINESTREBATTITO_ECG{i-iniziopicchi+1} = finestrabattito_ECG;
%     FINESTREBATTITO_SCG{i-iniziopicchi+1} = finestrabattito_SCG;
%     clearvars finestrabattito_ECG finestrabattito_SCG

    picchi(i-iniziopicchi+1,:) = n_picchi;
   
 end 
     name = erase(name,"ECG-FILT-");
    save(['C:\Users\feder\Desktop\Tesi\Data\Windows SCG' 'Windows SCG-' name],'picchi','FINESTREBATTITO_SCG','FINESTREBATTITO_ECG')
end 