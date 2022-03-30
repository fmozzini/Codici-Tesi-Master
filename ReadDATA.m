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
N = length(list)
list(N-1) = []; 
list(N-1) = [];
N = length(list); %N = numero acquisizioni = 8 da 24 ore
%%

for i = 1:N %Tutti i soggetti.    
    FOLDER = fullfile(list(i).folder, list(i).name);
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
    Ecg = [ Ecg Time_Ecg ];
    
    namefile = list(i).name

    save(['C:\Users\feder\Desktop\Tesi\Data\24h Signals\' namefile '.mat'],'Acc','Rot','Ecg')
    save(['C:\Users\feder\Desktop\Tesi\Data\ECG\' 'ECG-' namefile '.mat'],'Ecg')
    
end 

% %%
% % Ora capisco come fare il filtraggio
% % Z ACCELERATION
% z_filt = filtHP(Acc.Z_Acc,fs_Acc);
% 
% % figure()
% % subplot(211), plot(Acc.Z_Acc(1:640)), hold on, plot(z_filt(1:640))
% % subplot(212), plot(Ecg.Values(1:10240))
% % sgtitle('Z Acceleration - filterHP')
% figure()
% subplot(211), plot(Acc.Z_Acc(6*230400:6*230400+320)), hold on, plot(z_filt(6*230400:6*230400+320))
% subplot(212), plot(Ecg.Values(6*3686400:6*3686400+5120))
% sgtitle('Z Acceleration - filtroHP')
% 
% % X ACCELERATION
% x_filt = filtroR(Acc.X_Acc,fs_Acc);
% 
% % figure()
% % subplot(211), plot(Acc.X_Acc(1:640)), hold on, plot(x_filt(1:640))
% % subplot(212), plot(Ecg.Values(1:10240))
% % sgtitle('X Acceleration - filtroR')
% figure()
% subplot(211), plot(Acc.X_Acc(6*230400:6*230400+320)), hold on, plot(x_filt(6*230400:6*230400+320))
% subplot(212), plot(Ecg.Values(6*3686400:6*3686400+5120))
% sgtitle('X Acceleration - filtroR')
% 
% % Y ACCELERATION
% y_filt = filtroR(Acc.Y_Acc,fs_Acc);
% 
% % 
% % figure()
% % subplot(211), plot(Acc.Y_Acc(1:640)), hold on, plot(y_filt(1:640))
% % subplot(212), plot(Ecg.Values(1:10240))
% % sgtitle('Y Acceleration - filtroR')
% figure()
% subplot(211), plot(Acc.Y_Acc(6*230400:6*230400+320)), hold on, plot(y_filt(6*230400:6*230400+320))
% subplot(212), plot(Ecg.Values(6*3686400:6*3686400+5120))
% sgtitle('Y Acceleration - filtroR')
% 
% % save('x_filt'); save('y_filt'); save('z_filt')

%%
% X - ROTATION 
% x_rotfilt = filtroR(Rot.X_Rot,fs_Rot);
% 
% % figure()
% % subplot(211), plot(Rot.X_Rot(1:640)), hold on, plot(x_rotfilt(1:640))
% % subplot(212), plot(Ecg.Values(1:10240))
% % sgtitle('X Rotation - filtroR')
% 
% figure()
% subplot(211), plot(Rot.X_Rot(6*230400:6*230400+320)), hold on, plot(x_rotfilt(6*230400:6*230400+320))
% subplot(212), plot(Ecg.Values(6*3686400:6*3686400+5120))
% sgtitle('X Rotation - filtroR')
% 
% % Y ROTATION
% y_rotfilt = filtroR(Rot.Y_Rot,fs_Rot);
% 
% % figure()
% % subplot(211), plot(Rot.Y_Rot(1:640)), hold on, plot(y_rotfilt(1:640))
% % subplot(212), plot(Ecg.Values(1:10240))
% % sgtitle('Y Rotation - filtroR')
% 
% figure()
% subplot(211), plot(Rot.Y_Rot(6*230400:6*230400+320)), hold on, plot(y_rotfilt(6*230400:6*230400+320))
% subplot(212), plot(Ecg.Values(6*3686400:6*3686400+5120))
% sgtitle('Y Rotation - filtroR')
% 
% % Z ROTATION
% z_rotfilt = filtHP(Rot.Z_Rot,fs_Rot);
% 
% % figure()
% % subplot(211), plot(Rot.Z_Rot(1:640)), hold on, plot(z_rotfilt(1:640))
% % subplot(212), plot(Ecg.Values(1:10240))
% % sgtitle('Z Rotation - filterHP')
% 
% figure()
% subplot(211), plot(Rot.Z_Rot(6*230400:6*230400+320)), hold on, plot(z_rotfilt(6*230400:6*230400+320))
% subplot(212), plot(Ecg.Values(6*3686400:6*3686400+5120))
% sgtitle('Z Rotation - filtroHP')

% %%
% % Acceleration
%      f1 = figure()
%      sgtitle('Acceleration x,y,z VS Filtered Acceleration')
%      subplot(321), plot(Acc.Time_Acc./3600, Acc.X_Acc,'r'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(x)'), axis tight
%      subplot(322), plot(Acc.Time_Acc./3600, x_filt,'r'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(x) FILTERED'), axis tight
%      subplot(323), plot(Acc.Time_Acc./3600, Acc.Y_Acc,'g'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(y)'), axis tight
%      subplot(324), plot(Acc.Time_Acc./3600, y_filt,'g'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(y) FILTERED'), axis tight
%      subplot(325), plot(Acc.Time_Acc./3600, Acc.Z_Acc,'b'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(z)'), axis tight
%      subplot(326), plot(Acc.Time_Acc./3600, z_filt,'b'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(z) FILTERED'), axis tight
%      
%      % Rotation
%      f2 = figure()
%      sgtitle('Angular Rate x,y,z VS Filtered Angular Rate')
%      subplot(321), plot(Rot.Time_Rot./3600, Rot.X_Rot,'r'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(x)'), axis tight
%      subplot(322), plot(Rot.Time_Rot./3600, x_rotfilt,'r'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(x) FILTERED'), axis tight
%      subplot(323), plot(Rot.Time_Rot./3600, Rot.Y_Rot,'g'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(y)'), axis tight
%      subplot(324), plot(Rot.Time_Rot./3600, y_rotfilt,'g'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(y) FILTERED'), axis tight
%      subplot(325), plot(Rot.Time_Rot./3600, Rot.Z_Rot,'b'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(z)'), axis tight
%      subplot(326), plot(Rot.Time_Rot./3600, z_rotfilt,'b'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(z) FILTERED'), axis tight
%% 
% COMUNQUE TROPPO GRANDIIII! 
% figure()
% subplot(311), plot(x_filt_new(6*230400:6*230400+320))
% subplot(312), plot(y_filt_new(6*230400:6*230400+320))
% subplot(313), plot(z_filt_new(6*230400:6*230400+320))

% cinetical = KE(55,x_filt,y_filt,z_filt,Acc.Time_Acc);
% %cineticar = KE(55,x_rotfilt,y_rotfilt,z_rotfilt,Rot.Time_Rot);
% 
%  % all'ora 64*60*60 = 230400
%  % all'ora 1024*60*60 = 3686400
% % lineare 
% figure()
% subplot(211), plot(cinetical(6*230400:6*230400+320))
% subplot(212), plot(Ecg.Values(6*3686400:6*3686400+5120)),
% sgtitle('Kinetic Energy Lineare')
% 
% cineticar = KEr(53,165,x_rotfilt,y_rotfilt,z_rotfilt,Rot.Time_Rot);
% 
% figure()
% subplot(311), plot(cinetical(6*230400:6*230400+320),'r'),title('Kinetic Energy Linear')
% subplot(312), plot(cineticar(6*230400:6*230400+320),'g'),title('Kinetic Energy Rotational')
% subplot(313), plot(Ecg.Values(6*3686400:6*3686400+5120)),title('Ecg')
% 
