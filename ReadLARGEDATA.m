% Filtered signals 
clc
close all
clear all

%% Read data 
prompt = {'Seleziona cartella:'};
name = '';
numlines = 1;
defaultanswer = {'Movisens acquisitions'};
% cambio il nome per ogni paziente -> ho 8 cartelle, in questo modo ho sempre
% che il primo elemento è acc ed il secondo è rot
answer = inputdlg(prompt, name, numlines, defaultanswer); %crea finestra di dialogo per selezione dei dati
ris = fullfile('C:\Users\feder\Desktop\Tesi\',char(answer));

list = dir(ris); %lista dei file contenuti nella cartella con relative path=ris
list(1) = []; %cancella la prima riga --> the current folder: '.'
list(1) = [];%cancella di nuovo la prima riga --> the 'up a level' folder: '..'
N = length(list);
list(N-1) = []; 
list(N-1) = [];
N = length(list)
%%
for i = N:N
    FOLDER = fullfile(list(i).folder, list(i).name);
    file = dir(FOLDER);
    file(1) = [];
    file(1) = [];
    Acc = openLargeDatafiles(fullfile(file(1).folder,file(1).name));
    fs_Acc = 64;
    Rot = openLargeDatafiles(fullfile(file(2).folder,file(2).name));
    fs_Rot = 64;
    Ecg = openLargeDatafiles(fullfile(file(5).folder,file(5).name));
    fs_Ecg = 1024;

    % Prendo solo le 24 ore che mi interessano, devo cercare il punto di
    % partenza
    GiornoSCG = 24*3600*fs_Acc;
    GiornoECG = 24*3600*fs_Ecg;
    inizio_campioni_SCG = 4746800;
    inizio_campioni_ECG = (inizio_campioni_SCG/fs_Acc)*fs_Ecg;
    Acc = Acc(inizio_campioni_SCG:inizio_campioni_SCG+GiornoSCG,:);
    Rot = Rot(inizio_campioni_SCG:inizio_campioni_SCG+GiornoSCG,:);
    Ecg = Ecg(inizio_campioni_ECG:inizio_campioni_ECG+GiornoECG,:);

    Acc = table(Acc(:,1),Acc(:,2),Acc(:,3));
    Acc.Properties.VariableNames = {'X_Acc' 'Y_Acc' 'Z_Acc'};
    Rot = table(Rot(:,1),Rot(:,2),Rot(:,3));
    Rot.Properties.VariableNames = {'X_Rot' 'Y_Rot' 'Z_Rot'};
    Ecg = table(Ecg);
    Ecg.Properties.VariableNames = {'Values'};

    Time_Acc = ((0:numel(Acc.X_Acc)-1)/fs_Acc)';
    Time_Acc = table(Time_Acc); %seconds
    Acc = [ Acc Time_Acc ]; %Acc.X_Acc Acc.Y_Acc Acc.Z_Acc Time_Acc[s]

    Time_Rot = ((0:numel(Rot.X_Rot)-1)/fs_Rot)';
    Time_Rot = table(Time_Rot); %seconds
    Rot = [ Rot Time_Rot ]; %Rot.X_Rot Rot.Y_Rot Rot.Z_Rot Time_Rot[s]
    
    Time_Ecg = ((0:numel(Ecg.Values)-1)/fs_Ecg)';
    Time_Ecg = table(Time_Ecg); %seconds
    Ecg = [ Ecg Time_Ecg ];

    namefile = list(i).name
    save(['C:\Users\feder\Desktop\Tesi\Data\24h Signals\' namefile '.mat'],'Acc','Rot','Ecg')
    save(['C:\Users\feder\Desktop\Tesi\Data\ECG\' 'ECG-' namefile '.mat'],'Ecg')
end 
figure()
subplot(211),plot(Acc(:,3))
subplot(212),plot(Acc1(:,3))

figure()
plot(Acc.Z_Acc)
