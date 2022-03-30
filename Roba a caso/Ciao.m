clear all
clc
close all

%%
prompt={'Seleziona cartella:'}; 
name='';
numlines=1;
defaultanswer={'2021-06-12 09.38.41 MAIOLO'};

answer=inputdlg(prompt,name,numlines,defaultanswer); %crea finestra di dialogo per selezione dei dati

%ris=char(answer);

ris = fullfile('C:\Users\feder\Desktop\Tesi\Movisens acquisitions',char(answer));
%ris = 'C:\Users\feder\Desktop\Tesi\Movisens acquisitions\2021-06-12 09.38.41 MAIOLO';
list = dir(ris); %lista dei file contenuti nella cartella con relative path=ris
list(1) = []; %cancella la prima riga --> the current folder: '.' 
list(1) = []; %cancella di nuovo la prima riga --> the 'up a level' folder: '..'

N = length(list);
% Loro avevano 90 acquistizione (3 acquisizioni per 30 persone)
% N = length(list); %N = 8 ( noi abbiamo 1 acqusizione per soggetto) -> 8
% cartelle
LIMIT = 58/100*N; %TRAINING SET = 58% DEL DATASET PER OGNI DIREZIONE X, Y, Z (58% di 90 è 53)
train_count = 0; %dimensione del training set
test_count = 0; %dimensione del test set
training_flag = 0; %flag che determina lo stato del training: 0=ongoing 1=end

%FOLDER=fullfile(list.folder,list.name);
FOLDER = ris;
file = dir(FOLDER);
file(1) = [];
file(1) = [];
%filename = strcat('/Users/feder/Desktop/Tesi/Movisens acquisitions/2021-06-12 09.38.41 MAIOLO/',...
%'angularrate.csv'); % INSERIRE PATH CORRETTO

%contatori per ogni direzione
x_count = 0;
y_count = 0;
z_count = 0;

    if training_flag == 0 % stato ongoing
%% Signals in 3 directions 
    Acc = readtable(fullfile(file(1).folder,file(1).name));
    Acc.Properties.VariableNames = {'X_Acc' 'Y_Acc' 'Z_Acc'};
    fs_Acc = 64;
    Rot = readtable(fullfile(file(2).folder,file(2).name));
    Rot.Properties.VariableNames = {'X_Rot' 'Y_Rot' 'Z_Rot'};
    fs_Rot = 64;

    Time_Acc = ((0:numel(Acc.X_Acc)-1)/fs_Acc)';
    Time_Acc = table(Time_Acc); %seconds
    Acc = [ Acc Time_Acc ]; %Acc.X_Acc Acc.Y_Acc Acc.Z_Acc Time_Acc[s]

    Time_Rot = ((0:numel(Rot.X_Rot)-1)/fs_Rot)';
    Time_Rot = table(Time_Rot); %seconds
    Rot = [ Rot Time_Rot ]; %Rot.X_Rot Rot.Y_Rot Rot.Z_Rot Time_Rot[s]
    %[Timestamp_Acc, X_Acc , Y_Acc , Z_Acc] = import_acc_csv(fullfile(file(1).folder,file(1).name)); % se il file è .xlxl usa funzione import_acc e import_rot
    %[Timestamp_Rot, X_Rot , Y_Rot , Z_Rot ] = import_rot_csv(fullfile(file(2).folder,file(2).name));   

    LIMIT = 58/100*N; %TRAINING SET = 58% DEL DATASET PER OGNI DIREZIONE X, Y, Z (58% di 90 è 53)
    train_count = 0; %dimensione del training set
    test_count = 0; %dimensione del test set
    training_flag = 0; %flag che determina lo stato del training: 0=ongoing 1=end

    %% Signal Portion Identification        
    v=4; %default threshold is 4*std
    [X_peak,X_loc] = findpeaks(x,'MinPeakHeight',v*std(x));
    [Y_peak,Y_loc] = findpeaks(y,'MinPeakHeight',v*std(y));
    [Z_peak,Z_loc] = findpeaks(z,'MinPeakHeight',v*std(z));

    loc = {X_loc Y_loc Z_loc};
    
    %f1=figure('Name',strcat('Record n.',num2str(i)), 'NumberTitle', 'off') , %% Accelerometer
    f1 = figure()
    f1.WindowState = 'maximized';
    plot(x),hold on, plot(X_loc,x(X_loc),'*r');