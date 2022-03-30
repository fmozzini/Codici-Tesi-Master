%% New code: checkbox 
% ANALIZZI PAZIENTE PER VOLTA ED ORA PER VOLTA (COSI' NON HO PROBLEMI NEL
% SALVATAGGIO). PLOTTO TUTTI E 6 I SEGNALI INSIEME CON BOX VICINO PER
% RISPOSTA (POSSO METTERE CHE DI DEFAULT IL SEGNALE E' PULITO), IN BASSO
% METTO PULSANTE PER DIRE CHE SONO TUTTI PULITI. IN QUESTO MODO POSSO
% PENSARE COSA FARE CON I MINUTI... MAGARI PRIMA ANALISI 
% CONSIDERO 3 LABEL: RUMORE, TUTTI E 10 I BATTITI, NON TUTTI 10 I BATTITI
% (QUANTI NE CONSIDERO) ???  BOHHHHHHHHHHHH
% Chiedere a Sara per energia cinetica ( potrebbe servire ). 

% CONTROLLARE IL SEGNALE FILTRATO SUI 10 SECONDI O QUALCOSA DEL GENERE! 

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
   
%    %% Plot of accelerometers ( no filtered)
%    f1 = figure()
%    sgtitle('Acceleration x,y,z SECONDS')
%    subplot(311), plot(Acc.Time_Acc, Acc.X_Acc,'r'),xlabel('Time[s]'), ylabel('Acc(x)[m/s^2]'),legend('acc(x)'),axis tight
%    subplot(312), plot(Acc.Time_Acc, Acc.Y_Acc,'g'),xlabel('Time[s]'), ylabel('Acc(y)[m/s^2]'),legend('acc(y)'),axis tight
%    subplot(313), plot(Acc.Time_Acc, Acc.Z_Acc,'b'),xlabel('Time[s]'), ylabel('Acc(z)[m/s^2]'),legend('acc(z)'),axis tight
%    
%    f2 = figure()
%    sgtitle('Acceleration x,y,z HOURS ')
%    subplot(311), plot(Acc.Time_Acc./3600, Acc.X_Acc,'r'),xlabel('Time[h]'), ylabel('Acc(x)[m/s^2]'),legend('acc(x)'),axis tight
%    subplot(312), plot(Acc.Time_Acc./3600, Acc.Y_Acc,'g'),xlabel('Time[h]'), ylabel('Acc(y)[m/s^2]'),legend('acc(y)'),axis tight
%    subplot(313), plot(Acc.Time_Acc./3600, Acc.Z_Acc,'b'),xlabel('Time[h]'), ylabel('Acc(z)[m/s^2]'),legend('acc(z)'),axis tight
%    
     %% Filtering with Butterworth bandpass filter (10-13 Hz) ACCELERATION
     % %% Butterworth bandpass filter of order 6 (0.5-25 Hz) --> X ACC
     n = 4;
     ftype = 'bandpass';
     Wn = [5 25]/(fs_Acc/2);
     [b,a] = butter(n,Wn,ftype);
     %fcuthigh=0.01; %KHz
     %fcutlow=0.013; %KHz
     %[b,a]=butter(6,[fcutlow,fcuthigh],'bandpass');
     Acc_xfilt = filter(b,a,Acc.X_Acc);
     Acc_yfilt = filter(b,a,Acc.Y_Acc);
     Acc_zfilt = filter(b,a,Acc.Z_Acc);
     DATI_ACC{i,1} = Acc_xfilt;
     DATI_ACC{i,2} = Acc_yfilt;
     DATI_ACC{i,3} = Acc_zfilt;
     
     %% Filtering with Butterworth bandpass filter (10-13 Hz) GYROSCOPE
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
     
     %% Comparison of filtered signal and normal signal
     % Acceleration
     f1 = figure()
     sgtitle('Acceleration x,y,z VS Filtered Acceleration')
     subplot(321), plot(Acc.Time_Acc./3600, Acc.X_Acc,'r'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(x)'), axis tight
     subplot(322), plot(Acc.Time_Acc./3600, Acc_xfilt,'r'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(x) FILTERED'), axis tight
     subplot(323), plot(Acc.Time_Acc./3600, Acc.Y_Acc,'g'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(y)'), axis tight
     subplot(324), plot(Acc.Time_Acc./3600, Acc_yfilt,'g'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(y) FILTERED'), axis tight
     subplot(325), plot(Acc.Time_Acc./3600, Acc.Z_Acc,'b'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(z)'), axis tight
     subplot(326), plot(Acc.Time_Acc./3600, Acc_zfilt,'b'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(z) FILTERED'), axis tight
     
     % Rotation
     f2 = figure()
     sgtitle('Rotation x,y,z VS Filtered Rotation')
     subplot(321), plot(Rot.Time_Rot./3600, Rot.X_Rot,'r'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(x)'), axis tight
     subplot(322), plot(Rot.Time_Rot./3600, Rot_xfilt,'r'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(x) FILTERED'), axis tight
     subplot(323), plot(Rot.Time_Rot./3600, Rot.Y_Rot,'g'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(y)'), axis tight
     subplot(324), plot(Rot.Time_Rot./3600, Rot_yfilt,'g'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(y) FILTERED'), axis tight
     subplot(325), plot(Rot.Time_Rot./3600, Rot.Z_Rot,'b'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(z)'), axis tight
     subplot(326), plot(Rot.Time_Rot./3600, Rot_zfilt,'b'), xlabel('Time[h]'), ylabel('Rot[rad/s]'),legend('rot(z) FILTERED'), axis tight
     
     % FAI ANALISI ANCHE SU 10 SEC
     
     % ANALISI DI UN'ORA
     samplesh = 64*3600;
     L = size(Acc.X_Acc,1);
     SecH = 3600;
     frlen = 10; %10 secondi
     H = 1; % Hour that we are analyzing
     frame_number = 0;
     
      f3 = figure()
      Hour_n = 0; % QUALE ORA VUOI GUARDARE?
      for samples = 1: samplesh:((1*samplesh)-1)
      subplot(6,2,1), plot((Acc.Time_Acc)./3600,DATI_ACC{i,1}), title('Acclearion(x) 24 hour'), xlabel('Time[h]'), ylabel('[m/s^2]')    
      subplot(6,2,2), plot((Acc.Time_Acc(samples:((samples+samplesh)-1))./3600), DATI_ACC{i,1}(samples:((samples+samplesh)-1))), title(['Accelarion(x) Hour n°' num2str(Hour_n+1)]), xlabel('Time[h]'), ylabel('[m/s^2]')
      subplot(6,2,3), plot((Acc.Time_Acc)./3600,DATI_ACC{i,2}), title('Acclearion(y) 24 hour'), xlabel('Time[h]'), ylabel('[m/s^2]')    
      subplot(6,2,4), plot((Acc.Time_Acc(samples:((samples+samplesh)-1))./3600), DATI_ACC{i,2}(samples:((samples+samplesh)-1))), title(['Accelarion(y) Hour n°' num2str(Hour_n+1)]), xlabel('Time[h]'), ylabel('[m/s^2]')
      subplot(6,2,5), plot((Acc.Time_Acc)./3600,DATI_ACC{i,3}), title('Acclearion(z) 24 hour'), xlabel('Time[h]'), ylabel('[m/s^2]')    
      subplot(6,2,6), plot((Acc.Time_Acc(samples:((samples+samplesh)-1))./3600), DATI_ACC{i,3}(samples:((samples+samplesh)-1))), title(['Accelarion(z) Hour n°' num2str(Hour_n+1)]), xlabel('Time[h]'), ylabel('[m/s^2]')
      subplot(6,2,7), plot((Rot.Time_Rot)./3600,DATI_ROT{i,1}), title('Rotation(x) 24 hour'), xlabel('Time[h]'), ylabel('[rad/s]')    
      subplot(6,2,8), plot((Rot.Time_Rot(samples:((samples+samplesh)-1))./3600), DATI_ROT{i,1}(samples:((samples+samplesh)-1))), title(['Rotation(x) Hour n°' num2str(Hour_n+1)]), xlabel('Time[h]'), ylabel('[rad/s]') 
      subplot(6,2,9), plot((Rot.Time_Rot)./3600,DATI_ROT{i,2}), title('Rotation(y) 24 hour'), xlabel('Time[h]'), ylabel('[rad/s]')    
      subplot(6,2,10), plot((Rot.Time_Rot(samples:((samples+samplesh)-1))./3600), DATI_ROT{i,2}(samples:((samples+samplesh)-1))), title(['Rotation(y) Hour n°' num2str(Hour_n+1)]), xlabel('Time[h]'), ylabel('[rad/s]') 
      subplot(6,2,11), plot((Rot.Time_Rot)./3600,DATI_ROT{i,3}), title('Rotation(z) 24 hour'), xlabel('Time[h]'), ylabel('[rad/s]')    
      subplot(6,2,12), plot((Rot.Time_Rot(samples:((samples+samplesh)-1))./3600), DATI_ROT{i,3}(samples:((samples+samplesh)-1))), title(['Rotation(z) Hour n°' num2str(Hour_n+1)]), xlabel('Time[h]'), ylabel('[rad/s]') 
      Hour_n = Hour_n +1 ;
      end 
      
     f4 = figure()  
     for j = 1 : frlen : (SecH-frlen)  
          sgtitle(['Hour n°',num2str(Hour_n)])
          subplot(331),  plot((Acc.Time_Acc(j*fs_Acc:(j+(frlen-1))*fs_Acc)./3600),DATI_ACC{i,1}(j*fs_Acc:(j+(frlen-1))*fs_Acc)) , title(['Acc x Frame n°' num2str(frame_number+1)]), xlabel('Time[h]'), ylabel('[m/s^2]')
          subplot(332),  plot((Acc.Time_Acc(j*fs_Acc:(j+(frlen-1))*fs_Acc)./3600),DATI_ACC{i,2}(j*fs_Acc:(j+(frlen-1))*fs_Acc)) , title(['Acc y Frame n°' num2str(frame_number+1)]), xlabel('Time[h]'), ylabel('[m/s^2]')
          subplot(333),  plot((Acc.Time_Acc(j*fs_Acc:(j+(frlen-1))*fs_Acc)./3600),DATI_ACC{i,3}(j*fs_Acc:(j+(frlen-1))*fs_Acc)) , title(['Acc z Frame n°' num2str(frame_number+1)]), xlabel('Time[h]'), ylabel('[m/s^2]')
          subplot(334),  plot((Rot.Time_Rot(j*fs_Rot:(j+(frlen-1))*fs_Rot)./3600),DATI_ROT{i,1}(j*fs_Rot:(j+(frlen-1))*fs_Rot)) , title(['Rot x Frame n°' num2str(frame_number+1)]), xlabel('Time[h]'), ylabel('[rad/2]')
          subplot(335),  plot((Rot.Time_Rot(j*fs_Rot:(j+(frlen-1))*fs_Rot)./3600),DATI_ROT{i,2}(j*fs_Rot:(j+(frlen-1))*fs_Rot)) , title(['Rot y Frame n°' num2str(frame_number+1)]), xlabel('Time[h]'), ylabel('[rad/2]')
          subplot(336),  plot((Rot.Time_Rot(j*fs_Rot:(j+(frlen-1))*fs_Rot)./3600),DATI_ROT{i,3}(j*fs_Rot:(j+(frlen-1))*fs_Rot)) , title(['Rot z Frame n°' num2str(frame_number+1)]), xlabel('Time[h]'), ylabel('[rad/2]')
          subplot(337),  plot((Ecg.Time_Ecg(j*fs_Ecg:(j+(frlen-1))*fs_Ecg)./3600),Ecg.Values(j*fs_Ecg:(j+(frlen-1))*fs_Ecg)), title(['Ecg Frame n°' num2str(frame_number+1)]), xlabel('Time[h]')
          subplot(338),  plot((Ecg.Time_Ecg(j*fs_Ecg:(j+(frlen-1))*fs_Ecg)./3600),Ecg.Values(j*fs_Ecg:(j+(frlen-1))*fs_Ecg)), title(['Ecg Frame n°' num2str(frame_number+1)]), xlabel('Time[h]')
          subplot(339),  plot((Ecg.Time_Ecg(j*fs_Ecg:(j+(frlen-1))*fs_Ecg)./3600),Ecg.Values(j*fs_Ecg:(j+(frlen-1))*fs_Ecg)), title(['Ecg Frame n°' num2str(frame_number+1)]), xlabel('Time[h]')
          frame_number = frame_number +1;
          
         
          Title = ['Analyze the Hour n°',num2str(Hour_n),' Frame n°', num2str(frame_number)];

            %%%% SETTING DIALOG OPTIONS
            % Options.WindowStyle = 'modal';
            Options.Resize = 'on';
            Options.Interpreter = 'tex';
            Options.CancelButton = 'on';
            %Options.ApplyButton = 'on';
            Options.ButtonNames = {'Continue','Cancel'}; %<- default names, included here just for illustration
            Option.Dim = 3; % Horizontal dimension in fields

            Prompt = {};
            Formats = {};
            DefAns = struct([]);

%             Prompt(1,:) = {['Are these signals motion corrupted? '],[],[]};
%             Formats(1,1).type = 'text';
%             Formats(1,1).size = [-1 0];
%             Formats(1,1).span = [1 2]; % item is 1 field x 4 fields

%             Prompt(2,:) = {'Hour Number', 'Hour',[]};
%             Formats(2,1).type = 'edit';
%             Formats(2,1).format = 'integer';
%             Formats(2,1).size = 200; % automatically assign the height
%             DefAns(1).Hour = Hour_n;

            Prompt(1,:) = {'ACC X','AccX',[]};
            Formats(1,1).type = 'list';
            Formats(1,1).format = 'text';
            Formats(1,1).style = 'radiobutton';
            Formats(1,1).items = {'Rest' 'Motion' 'Not Sure'};
            % Formats(7,1).span = [2 1];  % item is 2 field x 1 fields
            DefAns(1).AccX = 'Motion';%3; % yen

            Prompt(2,:) = {'ACC Y','AccY',[]};
            Formats(1,2).type = 'list';
            Formats(1,2).format = 'text';
            Formats(1,2).style = 'radiobutton';
            Formats(1,2).items = {'Rest' 'Motion' 'Not Sure'};
            % Formats(7,1).span = [2 1];  % item is 2 field x 1 fields
            DefAns.AccY = 'Motion';%3; % yen

            Prompt(3,:) = {'ACC Z','AccZ',[]};
            Formats(1,3).type = 'list';
            Formats(1,3).format = 'text';
            Formats(1,3).style = 'radiobutton';
            Formats(1,3).items = {'Rest' 'Motion' 'Not Sure'};
            % Formats(7,1).span = [2 1];  % item is 2 field x 1 fields
            DefAns.AccZ = 'Motion';%3; % yen

            Prompt(4,:) = {'ROT X','RotX',[]};
            Formats(2,1).type = 'list';
            Formats(2,1).format = 'text';
            Formats(2,1).style = 'radiobutton';
            Formats(2,1).items = {'Rest' 'Motion' 'Not Sure'};
            % Formats(7,1).span = [2 1];  % item is 2 field x 1 fields
            DefAns.RotX = 'Motion';%3; % yen

            Prompt(5,:) = {'ROT Y','RotY',[]};
            Formats(2,2).type = 'list';
            Formats(2,2).format = 'text';
            Formats(2,2).style = 'radiobutton';
            Formats(2,2).items = {'Rest' 'Motion' 'Not Sure'};
            % Formats(7,1).span = [2 1];  % item is 2 field x 1 fields
            DefAns.RotY = 'Motion';%3; % yen

            Prompt(6,:) = {'ROT Z','RotZ',[]};
            Formats(2,3).type = 'list';
            Formats(2,3).format = 'text';
            Formats(2,3).style = 'radiobutton';
            Formats(2,3).items = {'Rest' 'Motion' 'Not Sure' };
            % Formats(7,1).span = [2 1];  % item is 2 field x 1 fields
            DefAns.RotZ = 'Motion';%3; % yen
            [Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options)
            

            if (strcmp(Answer.AccX,'Motion'))
                choice{Hour_n,1}(j*fs_Acc:((j+frlen)*fs_Acc),:) = 1;
            elseif (strcmp(Answer.AccX,'Not Sure'))
                choice{Hour_n,1}(j*fs_Acc:((j+frlen)*fs_Acc),:) = 2;
            else
                choice{Hour_n,1}(j*fs_Acc:((j+frlen)*fs_Acc),:) = 0;
            end
            
            if (strcmp(Answer.AccY,'Motion'))
                choice{Hour_n,2}(j*fs_Acc:((j+frlen)*fs_Acc),:) = 1;
            elseif (strcmp(Answer.AccY,'Not Sure'))
                choice{Hour_n,2}(j*fs_Acc:((j+frlen)*fs_Acc),:) = 2;
            else
                choice{Hour_n,2}(j*fs_Acc:((j+frlen)*fs_Acc),:) = 0;
            end
            
            if (strcmp(Answer.AccZ,'Motion'))
                choice{Hour_n,3}(j*fs_Acc:((j+frlen)*fs_Acc),:) = 1;
            elseif (strcmp(Answer.AccZ,'Not Sure'))
                choice{Hour_n,3}(j*fs_Acc:((j+frlen)*fs_Acc),:) = 2;
            else
                choice{Hour_n,3}(j*fs_Acc:((j+frlen)*fs_Acc),:) = 0;
            end
            
            if (strcmp(Answer.RotX,'Motion'))
                choice{Hour_n,4}(j*fs_Rot:((j+frlen)*fs_Rot),:) = 1;
            elseif (strcmp(Answer.RotX,'Not Sure'))
                choice{Hour_n,4}(j*fs_Rot:((j+frlen)*fs_Rot),:) = 2;
            else
                choice{Hour_n,4}(j*fs_Rot:((j+frlen)*fs_Rot),:) = 0;
            end
            
            if (strcmp(Answer.RotY,'Motion'))
                choice{Hour_n,5}(j*fs_Rot:((j+frlen)*fs_Rot),:) = 1;
            elseif (strcmp(Answer.RotY,'Not Sure'))
                choice{Hour_n,5}(j*fs_Rot:((j+frlen)*fs_Rot),:) = 2;
            else
                choice{Hour_n,5}(j*fs_Rot:((j+frlen)*fs_Rot),:) = 0;
            end
            
            if (strcmp(Answer.RotZ,'Motion'))
                choice{Hour_n,6}(j*fs_Rot:((j+frlen)*fs_Rot),:) = 1;
            elseif (strcmp(Answer.RotZ,'Not Sure'))
                choice{Hour_n,6}(j*fs_Rot:((j+frlen)*fs_Rot),:) = 2;
            else
                choice{Hour_n,6}(j*fs_Rot:((j+frlen)*fs_Rot),:) = 0;
            end


     end 
           
     % SALVARE I FILE CHOICE -> FAI FILE MAT E ASSEGNA NOME CON INIZIO PERCORSO_ORA   
         
end 




figure()
plot(Ecg.Values(1:10240))
figure()
plot(Rot.Z_Rot(1:640))

% AGGIUSTA FILTRAGGIO ROT
% CALCOLA ENERGIA CINETICA PER SARAH PER ROT E ACC
% FACCIO LABEL CON CON ECG, ENERGIA CINETICA SU FINESTRE DI 10 S NON
% SOVRAPPOSTE. GUARDO 0-1 SOLO PER CI SONO TUTTI O NON CI SONO! 
% POTREBBE AVER SENSO METTERE UNA SOGLIA E FARE UN PRIMO RICONOSCIMENTO (O
% NO... BOHHH??? ) 