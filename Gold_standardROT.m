%% Gold standard - Acc
clc
close all
clear all

%% GOLD STANDARD FOR ALL DATASET
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
list(9) = []; 
list(9) = [];
N = length(list); %N = numero acquisizioni = 8 da 24 ore
%%
for i = 1:1 %Tutti i soggetti.    
    FOLDER = fullfile(list(1).folder, list(1).name);
    file = dir(FOLDER);
    file(1) = [];
    file(1) = [];
    
    Acc = readtable(fullfile(file(1).folder,file(1).name));
    Acc.Properties.VariableNames = {'X_Acc' 'Y_Acc' 'Z_Acc'};
    fs_Acc = 64;
    Rot = readtable(fullfile(file(2).folder,file(2).name));
    Rot.Properties.VariableNames = {'X_Rot' 'Y_Rot' 'Z_Rot'};
    fs_Rot = 64;
    Ecg = readtable(fullfile(file(5).folder,file(5).name));
    Ecg.Properties.VariableNames = {'Values'};
    fs_Ecg = 1024;

    Time_Acc = ((0:numel(Acc.X_Acc)-1)/fs_Acc)';
    Time_Acc = table(Time_Acc); %seconds
    Acc = [ Acc Time_Acc ]; %Acc.X_Acc Acc.Y_Acc Acc.Z_Acc Time_Acc[s]

    Time_Rot = ((0:numel(Rot.X_Rot)-1)/fs_Rot)';
    Time_Rot = table(Time_Rot); %seconds
    Rot = [ Rot Time_Rot ]; %Rot.X_Rot Rot.Y_Rot Rot.Z_Rot Time_Rot[s]
    
    Time_Ecg = ((0:numel(Ecg.Values)-1)/fs_Ecg)';
    Time_Ecg = table(Time_Ecg); %seconds
    Ecg = [ Ecg Time_Ecg ]; %Ecg.Values Time_Ecg[s]
   % sono stati importati i dati della misurazione
   % la frequenza di campionamento è costante, non devo fare un resampling
   
   %% Plot of rotations ( no filtered)
   f1 = figure()
   sgtitle('Rotation x,y,z SECONDS')
   subplot(311), plot(Rot.Time_Rot, Rot.X_Rot,'r'),xlabel('Time[s]'), ylabel('Rot(x)[rad/s]'),legend('rot(x)'),axis tight
   subplot(312), plot(Rot.Time_Rot, Rot.Y_Rot,'g'),xlabel('Time[s]'), ylabel('Rot(y)[rad/s]'),legend('rot(y)'),axis tight
   subplot(313), plot(Rot.Time_Rot, Rot.Z_Rot,'b'),xlabel('Time[s]'), ylabel('Rot(z)[rad/s]'),legend('rot(z)'),axis tight
   
   f2 = figure()
   sgtitle('Rotation x,y,z HOURS ')
   subplot(311), plot(Rot.Time_Rot./3600, Rot.X_Rot,'r'),xlabel('Time[h]'), ylabel('Rot(x)[rad/s]'),legend('rot(x)'),axis tight
   subplot(312), plot(Rot.Time_Rot./3600, Rot.Y_Rot,'g'),xlabel('Time[h]'), ylabel('Rot(y)[rad/s]'),legend('rot(y)'),axis tight
   subplot(313), plot(Rot.Time_Rot./3600, Rot.Z_Rot,'b'),xlabel('Time[h]'), ylabel('Rot(z)[rad/s]'),legend('rot(z)'),axis tight
   
   
     %% Filtering with Butterworth bandpass filter (10-13 Hz)
     % ORA FACCIO SOLO PER GIROSCOPIO 
     % %% Butterworth bandpass filter of order 6 (0.5-25 Hz) --> X ACC
     n = 3;
     ftype = 'bandpass';
     Wn = [10 13]/(fs_Rot/2);
     [b,a] = butter(n,Wn,ftype);
     %fcuthigh=0.01; %KHz
     %fcutlow=0.013; %KHz
     %[b,a]=butter(6,[fcutlow,fcuthigh],'bandpass');
     Rot_xfilt = filter(b,a,Rot.X_Rot);
     Rot_yfilt = filter(b,a,Rot.Y_Rot);
     Rot_zfilt = filter(b,a,Rot.Z_Rot);
     DATI_ROT{i,1} = Rot_xfilt;
     DATI_ROT{i,2} = Rot_yfilt;
     DATI_ROT{i,3} = Rot_zfilt;
%% Comparison of filtered signal and normal signal (rotation )
     f3 = figure()
     sgtitle('Rotation x,y,z VS Filtered Rotation')
     subplot(321), plot(Rot.Time_Rot./3600, Rot.X_Rot,'r'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(x)'), axis tight
     subplot(322), plot(Rot.Time_Rot./3600, Rot_xfilt,'r'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(x) FILTERED'), axis tight
     subplot(323), plot(Rot.Time_Rot./3600, Rot.Y_Rot,'g'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(y)'), axis tight
     subplot(324), plot(Rot.Time_Rot./3600, Rot_yfilt,'g'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(y) FILTERED'), axis tight
     subplot(325), plot(Rot.Time_Rot./3600, Rot.Z_Rot,'b'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(z)'), axis tight
     subplot(326), plot(Rot.Time_Rot./3600, Rot_zfilt,'b'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(z) FILTERED'), axis tight
     
     %% Signal manipulation for visualization
%     DATI_plot = DATI_ROT;
%     for j=1:3 %canali x,y,z
%         for k=1:size(DATI_ROT{i,j},1) %tutte i valori delle 24ore
%             if DATI_ROT{i,j}(k,1) > 0.1
%                 DATI_plot{i,j}(k,1) = 0.1;
%             end
%             if DATI_ROT{i,j}(k,1) < -0.1
%                 DATI_plot{i,j}(k,1) = 0.1;
%             end
%         end
%     end
     
     
     %% Framing parameters
     fs = fs_Rot;
     frlen = round(1*fs); %64 (ogni sec ho 64 campioni)
     hop = frlen;
        
     %% Manual classification in motion or rest frames
     for j=1:3 %numero canali (x,y,z)
        frame_number = 0;
        figure,
        
        for k=1:frlen:((size(DATI_ROT{i,j},1)-frlen)) %tutti i campioni ogni 64
             subplot(311) , plot(DATI_ROT{i,j},'b')
             hold on
             plot(k:k+(frlen-1),DATI_ROT{i,j}(k : k+(frlen-1)),'r'), title('Total signal')
             subplot(312),  plot(k:k+(frlen-1),DATI_ROT{i,j}(k : k+(frlen-1))) , title(['Frame n°' num2str(frame_number+1)])
              subplot(313),  plot(k:k+(fs_Ecg-1),Ecg.Values(k:k+(fs_Ecg-1))), title(['Ecg Frame n°' num2str(frame_number+1)])
             ButtonName = questdlg ('Motion?',...
                 'SAVE','Yes');
             if (strcmp(ButtonName, 'Yes'))
                 choice{i,j}(k:k+frlen,:) = 1;
             else
                 choice{i,j}(k:k+frlen,:) = 0;
             end
             frame_number = frame_number + 1;   
         end
         close all       
     end
end

 % PROVIAMO A PLOTTARE INSIEME IL SEGNALE NORMALE E QUELLO CON SAMPLING SU
 % FINESTRE DI 10 SEC
% choice: lunghezza in sample, leggermente diversa da DATI_GYR 
%(es. 16491 vs 16502 nel primo soggetto)