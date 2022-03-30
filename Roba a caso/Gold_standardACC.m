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
   
   %% Plot of accelerometers ( no filtered)
   f1 = figure()
   sgtitle('Acceleration x,y,z SECONDS')
   subplot(311), plot(Acc.Time_Acc, Acc.X_Acc,'r'),xlabel('Time[s]'), ylabel('Acc(x)[m/s^2]'),legend('acc(x)'),axis tight
   subplot(312), plot(Acc.Time_Acc, Acc.Y_Acc,'g'),xlabel('Time[s]'), ylabel('Acc(y)[m/s^2]'),legend('acc(y)'),axis tight
   subplot(313), plot(Acc.Time_Acc, Acc.Z_Acc,'b'),xlabel('Time[s]'), ylabel('Acc(z)[m/s^2]'),legend('acc(z)'),axis tight
   
   f2 = figure()
   sgtitle('Acceleration x,y,z HOURS ')
   subplot(311), plot(Acc.Time_Acc./3600, Acc.X_Acc,'r'),xlabel('Time[h]'), ylabel('Acc(x)[m/s^2]'),legend('acc(x)'),axis tight
   subplot(312), plot(Acc.Time_Acc./3600, Acc.Y_Acc,'g'),xlabel('Time[h]'), ylabel('Acc(y)[m/s^2]'),legend('acc(y)'),axis tight
   subplot(313), plot(Acc.Time_Acc./3600, Acc.Z_Acc,'b'),xlabel('Time[h]'), ylabel('Acc(z)[m/s^2]'),legend('acc(z)'),axis tight
   
     %% Filtering with Butterworth bandpass filter (10-13 Hz)
     % ORA FACCIO SOLO PER GIROSCOPIO 
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
     
  %% Comparison of filtered signal and normal signal (acceleration X)
     f3 = figure()
     sgtitle('Acceleration x,y,z VS Filtered Acceleration')
     subplot(321), plot(Acc.Time_Acc./3600, Acc.X_Acc,'r'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(x)'), axis tight
     subplot(322), plot(Acc.Time_Acc./3600, Acc_xfilt,'r'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(x) FILTERED'), axis tight
     subplot(323), plot(Acc.Time_Acc./3600, Acc.Y_Acc,'g'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(y)'), axis tight
     subplot(324), plot(Acc.Time_Acc./3600, Acc_yfilt,'g'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(y) FILTERED'), axis tight
     subplot(325), plot(Acc.Time_Acc./3600, Acc.Z_Acc,'b'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(z)'), axis tight
     subplot(326), plot(Acc.Time_Acc./3600, Acc_zfilt,'b'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),legend('acc(z) FILTERED'), axis tight
     
     for j = 1:3 %--> cambia tutto in DATA_ACQ
         samplesh = 64*3600;
         Hour_n = 0;
         loc = 2;
         f4 = figure;
         for samples = 1 : samplesh:((8*samplesh)-1)
             loc = loc+1;
             subplot(5,2,[1,2]), plot(Acc.Time_Acc./3600, DATI_ACC{i,j},'r'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),title('Acceleration(x) FILTERED'), axis tight
             subplot(5,2,loc)
             plot(samples:samples+(samplesh-1),DATI_ACC{i,j}(samples:samples+(samplesh-1))),xlabel('Samples'),ylabel('Acc(x)[m/s^2]'), title(['Hour n°' num2str(Hour_n+1)])
             Hour_n = Hour_n + 1;
         end   

         Hour_n = 8;
         loc = 2;
         f5 = figure;
         for samples = samplesh*8 : samplesh:((16*samplesh)-1)
             loc = loc+1;
             subplot(5,2,[1,2]), plot(Acc.Time_Acc./3600, DATI_ACC{i,j},'r'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),title('Acceleration(x) FILTERED'), axis tight
             subplot(5,2,loc)
             plot(samples:samples+(samplesh-1),DATI_ACC{i,j}(samples:samples+(samplesh-1))),xlabel('Samples'),ylabel('Acc(x)[m/s^2]'), title(['Hour n°' num2str(Hour_n+1)])
             Hour_n = Hour_n + 1;      
         end 

         Hour_n = 16;
         loc = 2;
         f6 = figure;
         for samples = samplesh*16 : samplesh:((24*samplesh)-1)
             loc = loc+1;
             subplot(5,2,[1,2]), plot(Acc.Time_Acc./3600, DATI_ACC{i,j},'r'), xlabel('Time[h]'), ylabel('Acc[m/s^2]'),title('Acceleration(x) FILTERED'), axis tight
             subplot(5,2,loc)
             plot(samples:samples+(samplesh-1),DATI_ACC{i,j}(samples:samples+(samplesh-1))),xlabel('Samples'),ylabel('Acc(x)[m/s^2]'), title(['Hour n°' num2str(Hour_n+1)])
             Hour_n = Hour_n + 1;
         end 
        
       ButtonName = questdlg ('Do you want to label as motion corrupted some hours?',...
                      'SAVE','Yes');
       if (strcmp(ButtonName, 'Yes'))
           prompt = {'Enter Motion corrupted hours (space-separated numbers):'};
           dlgtitle = 'Motion Corrupted Hours';
           dims = [1 50];
           answer1 = inputdlg(prompt, dlgtitle, dims);
           answer1 = answer1{1};
           HoursMotionCorrupted = str2num(answer1);
           SechAcc = 3600*64;
           for a = 1:(size(HoursMotionCorrupted,2))
               value = HoursMotionCorrupted(a);
               choice{i,j}((SechAcc*(value-1))+1:(SechAcc*value)-1,:) = 2;
           end 
           
%            ButtonName1 = questdlg('Do you want to analyise 10 min windows?',...
%                'SAVE','Yes');
%            if (strcmp(ButtonName1,'Yes'))
%                prompt1 = {'Which hours are you intrested in (space-separated numbers):'};
%                dlgtitle1 = 'Analyse hours';
%                dims1 = [1 50];
%                answer2 = inputdlg(prompt1, dlgtitle1, dims1);
%                answer2 = answer2{1};
%                HoursAnalyse = str2num(answer2);
%                for h = 1:(size(HoursAnalyse,2))
%                    h = HoursAnalyse 
           frame_number = 0;
           frlen = 10;
           Sech = 3600;
           for k=1:frlen:((size(DATI_ACC{i,j},1)-frlen))
               if (choice{i,j}(k*fs_Acc:((k+frlen)*fs_Acc),:) ~= 2)
                   subplot(311) , plot(Acc.Time_Acc./3600,DATI_ACC{i,j},'b')
                   hold on
                   plot((Acc.Time_Acc(k*fs_Acc:(k+(frlen-1))*fs_Acc)./3600),DATI_ACC{i,j}(k*fs_Acc:(k+(frlen-1))*fs_Acc),'r'), title('Total signal')
                   subplot(312),  plot((Acc.Time_Acc(k*fs_Acc:(k+(frlen-1))*fs_Acc)./3600),DATI_ACC{i,j}(k*fs_Acc:(k+(frlen-1))*fs_Acc)) , title(['Acc Frame n°' num2str(frame_number+1)])
                   subplot(313),  plot((Ecg.Time_Ecg(k*fs_Ecg:(k+(frlen-1))*fs_Ecg)./3600),Ecg.Values(k*fs_Ecg:(k+(frlen-1))*fs_Ecg)), title(['Ecg Frame n°' num2str(frame_number+1)])
                   ButtonName = questdlg ('Motion?',...
                    'SAVE','Yes');
                   if (strcmp(ButtonName, 'Yes'))
                       choice{i,j}(k*fs_Acc:((k+frlen)*fs_Acc),:) = 1;
                   else
                       choice{i,j}(k*fs_Acc:((k+frlen)*fs_Acc),:) = 0;
                   end
                   frame_number = frame_number + 1; 
                   close all
               end
           end 
            
       else
           frame_number = 0;
           frlen = 10;
           Sech = 3600;
           for k=1:frlen:((size(DATI_ACC{i,j},1)-frlen))
               subplot(311) , plot(Acc.Time_Acc./3600,DATI_ACC{i,j},'b')
               hold on
               plot((Acc.Time_Acc(k*fs_Acc:(k+(frlen-1))*fs_Acc)./3600),DATI_ACC{i,j}(k*fs_Acc:(k+(frlen-1))*fs_Acc),'r'), title('Total signal')
               subplot(312),  plot((Acc.Time_Acc(k*fs_Acc:(k+(frlen-1))*fs_Acc)./3600),DATI_ACC{i,j}(k*fs_Acc:(k+(frlen-1))*fs_Acc)) , title(['Acc Frame n°' num2str(frame_number+1)])
               subplot(313),  plot((Ecg.Time_Ecg(k*fs_Ecg:(k+(frlen-1))*fs_Ecg)./3600),Ecg.Values(k*fs_Ecg:(k+(frlen-1))*fs_Ecg)), title(['Ecg Frame n°' num2str(frame_number+1)])
               ButtonName = questdlg ('Motion?',...
                 'SAVE','Yes');
               if (strcmp(ButtonName, 'Yes'))
                   choice{i,j}(k*fs_Acc:((k+frlen)*fs_Acc),:) = 1;
               else
                   choice{i,j}(k*fs_Acc:((k+frlen)*fs_Acc),:) = 0;
               end
               frame_number = frame_number + 1;   
           end
           close all       
    end
  end 
end 

    

          
     %% Framing parameters
%      fs = fs_Rot;
%      %frlen = round(1*fs)*10; %64 (ogni sec ho 64 campioni -> finestra di 640 campioni 
%      frlen = 10;
%      hop = frlen;
%      samplesh = 64*3600;  
%      samplesm = 64*60;
%      
%      %% Manual classification - HOURS 
%      for j = 1:3 %numero di canali x,y,z
%          windowmin_number = 5;
%          figure,
%          frlen = 10;
%          Sech = 3600;
%          for w=1:frlen*60:(Sech-1)
%             subplot(311) , plot(Acc.Time_Acc./3600,DATI_ACC{i,j},'b')
%             hold on
%             plot((Acc.Time_Acc(w*fs_Acc:(w+(frlen-1))*fs_Acc)./3600),DATI_ACC{i,j}(w*fs_Acc:(w+(frlen-1))*fs_Acc),'r'), title('Total signal')
%             subplot(312),  plot((Acc.Time_Acc(w*fs_Acc:(w+(frlen-1))*fs_Acc)./3600),DATI_ACC{i,j}(w*fs_Acc:(w+(frlen-1))*fs_Acc)) , title(['Acc Window min°' num2str(windowmin_number*2)])
%             subplot(313),  plot((Ecg.Time_Ecg(w*fs_Ecg:(w+(frlen-1))*fs_Ecg)./3600),Ecg.Values(w*fs_Ecg:(w+(frlen-1))*fs_Ecg)), title(['Ecg Windown min°' num2str(windowmin_number*2)])
%          end 
%      
% %      
%      %% Manual classification in motion or rest frames
%      
%         
%         % visualizzo in tempo -> controllo finestra campioni 
%         %for k=1:frlen:((size(DATI_ACC{i,j},1)-frlen)) 
%         % in 10 secondi ho 640 campioni di SCG e 10240 di ECG 
%         % faccio frlen = 10
%         % k = contatore sui secondi -> Sech = 3600 -> 3600*4 = 14400
%         % LE FINESTRE SONO DI 9 SECONDI, NON DI 10 
%         frame_number = 0;
%         frlen = 10
%         Sech = 3600;
%         for k=Sech*4:frlen:((size(DATI_ACC{i,j},1)-frlen)) 
%              subplot(311) , plot(Acc.Time_Acc./3600,DATI_ACC{i,j},'b')
%              hold on
%              plot((Acc.Time_Acc(k*fs_Acc:(k+(frlen-1))*fs_Acc)./3600),DATI_ACC{i,j}(k*fs_Acc:(k+(frlen-1))*fs_Acc),'r'), title('Total signal')
%              subplot(312),  plot((Acc.Time_Acc(k*fs_Acc:(k+(frlen-1))*fs_Acc)./3600),DATI_ACC{i,j}(k*fs_Acc:(k+(frlen-1))*fs_Acc)) , title(['Acc Frame n°' num2str(frame_number+1)])
%              subplot(313),  plot((Ecg.Time_Ecg(k*fs_Ecg:(k+(frlen-1))*fs_Ecg)./3600),Ecg.Values(k*fs_Ecg:(k+(frlen-1))*fs_Ecg)), title(['Ecg Frame n°' num2str(frame_number+1)])
%              ButtonName = questdlg ('Motion?',...
%                  'SAVE','Yes');
%              if (strcmp(ButtonName, 'Yes'))
%                  choice{i,j}(k:k+frlen,:) = 1;
%              else
%                  choice{i,j}(k:k+frlen,:) = 0;
%              end
%              frame_number = frame_number + 1;   
%          end
%          close all       
%      end
%   
%      end 
% end 

% ANALIZZI PAZIENTE PER VOLTA ED ORA PER VOLTA (COSI' NON HO PROBLEMI NEL
% SALVATAGGIO). PLOTTO TUTTI E 6 I SEGNALI INSIEME CON BOX VICINO PER
% RISPOSTA (POSSO METTERE CHE DI DEFAULT IL SEGNALE E' PULITO), IN BASSO
% METTO PULSANTE PER DIRE CHE SONO TUTTI PULITI. IN QUESTO MODO POSSO
% PENSARE COSA FARE CON I MINUTI... MAGARI PRIMA ANALISI 
% CONSIDERO 3 LABEL: RUMORE, TUTTI E 10 I BATTITI, NON TUTTI 10 I BATTITI
% (QUANTI NE CONSIDERO) ???  BOHHHHHHHHHHHH
% Chiedere a Sara per energia cinetica ( potrebbe servire ). 

% CONTROLLARE IL SEGNALE FILTRATO SUI 10 SECONDI O QUALCOSA DEL GENERE! 