% Ricalcoliamo gli intervalli RR e guardiamo cosa succede in ECG quando RR
% è strano

clear all
close all
clc
%%
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

 addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Filtered ECG'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Data\Pan-Tompkins\
 addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Picchi SCG - Acc z'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Codes\

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

    fs_ECG = 1024;
    fs_SCG = 64;
    window_SCG = 0.2*fs_SCG;
    window_ECG = 0.2*fs_ECG;
    R = zeros(peaksFORwindow(end-1),2);

    for i = 1:peaksFORwindow(end-1)
        R(i,:) = [qrs_I(i) qrs_AMP(i)];
    end
   
    RR = diff(R(:,1)); %indici, non s o ms (divido per 1024)
    RR_sec = RR./1024; % lo plotto in sec sulle y, 84163
    RR_tempo = R(1:end-1,1);
%     RR_tempo = RR_tempo./1024;
    plot(RR_tempo,RR_sec)

%% CERCHIAMO UN MODO PER TROVARE GLI OUTLIERS -> elimino i picchi se sono diversi dal 30% della media tra valore precedente e successivo
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6931967/ 
    eliminare = 0;
        for r = 1:length(RR_sec)-2 
            media = (RR_sec(r)+RR_sec(r+2))/2;
            diff = media - RR_sec(r+1);
            realdiff = diff/media;
            if (abs(realdiff)) >= 0.3
                eliminare = eliminare + 1;
                ELIMINARE(eliminare,:) = [RR_sec(r+1) RR_tempo(r+1) r+1];
            end 
        end 
    RR_pp = [RR_sec RR_tempo];
    R_post_processing = R;
    for n = length(ELIMINARE):-1:1
    %     R_post_processing(find(R_post_processing ==  ELIMINARE(n,2)),:) = [];
        R_post_processing(ELIMINARE(n,3),:) = [];
        RR_pp(ELIMINARE(n,3),:) = [];
    end 
   % ne ho eliminati 480 ma rimane qualcosina (pochi pochi outliers)
   figure()
   plot(RR_pp(2:end-1,2),RR_pp(2:end-1,1)),xlabel('[s]'),ylabel('RR[s]'),title('RR interval POST PROCESSING'); 

%%
    outliers = isoutlier(RR_pp(:,1)); % individuo altri 6 outliers, volendo posso vedere cosa sono
    out = 0;
    for i = 1:length(outliers)
        if outliers(i) == 1
            out = out+1;
            RR_pp(i,3) = 1;
        end 
    end 
    righe_out = find(RR_pp(:,3)==1);
    valori_out = RR_pp(righe_out,2);
    ampiezze_out = RR_pp(righe_out,1);
    PUNTI = [valori_out ampiezze_out]; 

%     figure()
%     plot(RR_tempo,RR_sec); hold on; plot(PUNTI(:,1),PUNTI(:,2),'*r')
% 
    mediano = median(RR_pp(:,1));
    mediano = mediano*1024;
    ectopico = zeros(length(RR_pp),2);
    n_ect = 0;
    R_post_processing2 = R_post_processing;
    nuovo_picco = 0;
    deleted = 0;
    
%%
    for i = 2:length(PUNTI)
        succ = find(R==PUNTI(i)); % riga, non numero di campioni
        finestrabattito_ECG = ECG_filt(PUNTI(i,1)-window_ECG:R(succ+1)-window_ECG)';
        plot((PUNTI(i,1)-window_ECG:R(succ+1)-window_ECG),finestrabattito_ECG); hold on; xline(PUNTI(i,1)-window_ECG+mediano,'--r'),title(i)
        disp(succ)
        answer = questdlg('What do you whant to do?', ...
	    'Analysing errors ', ...
	    'Missed R peak','Ectopic beat','Delete','Missed R peak');

        switch answer
            case 'Missed R peak'
                 nuovo_picco = nuovo_picco + 1;
                 [BATTITIAMP, BATTITIPOS] = findpeaks(finestrabattito_ECG,'NPeaks',2,'SortStr','descend');
                 battiti = [BATTITIPOS' BATTITIAMP'];
                 battiti = sortrows(battiti,1);
                 figure()
                 R_new(nuovo_picco,:) = [PUNTI(i,1)-window_ECG+battiti(2,1) battiti(2,2)]; 
                 plot((PUNTI(i,1)-window_ECG:R(succ+1)-window_ECG),finestrabattito_ECG); hold on; plot(R_new(nuovo_picco,1),R_new(nuovo_picco,2),'*r'),title(i)

                 answer2 = questdlg('Are you sure?',...)
                     'Last check', ...
                     'Yes','No','Yes');
                 switch answer2
                     case 'Yes'
                     case 'No'
                         R_new(nuovo_picco,:) = [];
                         nuovo_picco = nuovo_picco - 1;
                         R_post_processing2(find(R_post_processing == PUNTI(i)),:) = [];
                         deleted = deleted + 1;
                         disp('Deleted')
                 end 

            case 'Ectopic beat'
                 n_ect = n_ect + 1;
                 ectopico(succ,1) = PUNTI(i);
                 ectopico(succ,2) = 1; % ATTENZIONE, NON E' IL NUMERO DI CAMPIONE MA LA RIGA 
                 % non devo aggiungere niente, questo battito lo ho già 

            case 'Delete'
%                   R_post_processing(find(R == PUNTI(i)),:) = []; % qui sbaglia qualcosa... pensa un attimo!!
                  deleted = deleted + 1;
                  R_deleted(deleted) = find(R_post_processing==PUNTI(i));
                  disp('Deleted')
                  % meglio se mi segno solo le righe che voglio eliminare,
                  % poi fuori faccio un ciclo for che parte dalla fine e
                  % man a mano cancello tutto! In questo modo non dovrei
                  % avere problemi! 
        end
       
        pause 
        close all
    end 
%%
        RR_pp2 = RR_pp;
        for d = deleted:-1:1 
            % in R_deleted ho le righe che voglio eliminare da R
            R_post_processing2(R_deleted(d),:) = []; 
            RR_pp2(R_deleted(d),:) = [];
        end 
        
        if exist('R_new','var')~=0
            R_post_processingfin = [R_post_processing2; R_new]; 
        end 
        R_post_processingfin = sortrows(R_post_processing2,1);
        % non vado ad aggiungere gli intervalli RR nuovi 

        % SALVO IL RISULTATO
%         name = erase(name,"ECG_FILT-")
%         save(['C:\Users\feder\Desktop\Tesi\Data\PostProc PT\' 'PostProc PT-' name],'R_post_processingfin','RR_pp2','peaksFORwindow')
 end 

        %% RIFACCIAMO L'ANALISI (ultima finale finale)
        
        figure()
        a=subplot(211),plot(RR_pp2(2:end,2),RR_pp2(2:end,1)),xlabel('[s]'),ylabel('RR[s]'),title('RR interval POST PROCESSING'); 
        b=subplot(212),boxplot(RR_pp2(2:end,1)),title('RR interval POST PROCESSING')

        figure()
        a = subplot(211),plot(RR_tempo(2:end,:),RR_sec(2:end,:)); hold on; plot(PUNTI(:,1),PUNTI(:,2),'*r'),title('RR before post_processing')
        b = subplot(212), plot(RR_pp2(2:end,2),RR_pp2(2:end,1)),xlabel('[s]'),ylabel('RR[s]'),title('RR interval POST PROCESSING'); 
        linkaxes([a b],'x')