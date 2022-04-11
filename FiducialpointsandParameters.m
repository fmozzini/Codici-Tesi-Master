%% Define fiducial points and extract paramerers - SCG 
%% Qui voglio estrarre l'SCG solo se è presente anche l'ECG.
 clear all
 close all
 clc
%%
folderSCG = 'C:\Users\feder\Desktop\Tesi\Data\Picchi SCG - Acc z';
folderPT = 'C:\Users\feder\Desktop\Tesi\Data\Pan-Tompkins'; 
folderECG = 'C:\Users\feder\Desktop\Tesi\Data\Filtered ECG';
% folderECG = 'C:\Users\feder\Desktop\Tesi\Data\PostProc ECG';
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
%  addpath 'C:\Users\feder\Desktop\Tesi'\Data\''\'PostProc ECG'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Data\Pan-Tompkins\
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

    cinqueminuti = 5*60 %300 sec 
    cinqueminuti_ECG = cinqueminuti*fs_ECG; %307200
    cinqueminuti_SCG = cinqueminuti*fs_SCG; %19200
    n_finestre = cinqueminuti/30;
    iniziopicchi = peaksFORwindow(n_finestre);
    numero_picchi_totale = peaksFORwindow(end-1)-iniziopicchi;
 
% Fiducial Points
    MC = zeros(peaksFORwindow(end-1),2); % in realtà non mi serve
    RE = zeros(peaksFORwindow(end-1),2); % in realtà non mi serve
    AO = zeros(peaksFORwindow(end-1),2);
    AC = zeros(peaksFORwindow(end-1),2);
    IVC = zeros(peaksFORwindow(end-1),2); 
    R = zeros(peaksFORwindow(end-1),2);
% Parameters from SCG
    t_IVCAO = zeros(peaksFORwindow(end-1),1);
    t_IVCAC = zeros(peaksFORwindow(end-1),1);
    amp_IVCAO = zeros(peaksFORwindow(end-1),1);
    amp_IVCAC = zeros(peaksFORwindow(end-1),1);
    slope_IVCAO = zeros(peaksFORwindow(end-1),1);

    FINESTREBATTITO_ECG = zeros(peaksFORwindow(end-1),2000);
    FINESTREBATTITO_SCG = zeros(peaksFORwindow(end-1),100);
    % prova per trovare la lunghezza massima?
    maxlength_SCG = 0;
    maxlength_ECG = 0;
    count = 0;
    
    
 %%
%  for i = 1:peaksFORwindow(end-1) %scorro tutti i picchi ECG
    for i = iniziopicchi:peaksFORwindow(end-1)-1 %scorro tutti i picchi ECG

    qrs1 = (qrs_I(i)/1024)*64;
    qrs2 = (qrs_I(i+1)/1024)*64;
    finestra_dopo = qrs2-window_SCG %cambiato
    finestra_prima = qrs1-window_SCG
    [row]=find(POS_picchi_SCG<finestra_dopo & POS_picchi_SCG>finestra_prima)
    n_picchi = size(row,1)
    R(i,:) = [qrs_I(i) qrs_AMP(i)];

    if (n_picchi == 0) || (n_picchi >= 3)
        % tutti i parameters per questa finestra sono NaN
        t_IVCAO(i,1) = NaN;
        t_IVCAC(i,1) = NaN;
        amp_IVCAO(i,1) = NaN;
        amp_IVCAC(i,1) = NaN;
        slope_IVCAO(i,1) = NaN;

    else %% n_picchi == 1 or n_picchi == 2
%         count = count+1;
        picchi = POS_picchi_SCG(row)
        finestrabattito_ECG = ECG_filt(qrs_I(i)-window_ECG:qrs_I(i+1)-window_ECG)';
        finestrabattito_SCG = Acc_z(qrs1-window_SCG:qrs2-window_SCG)';

%         FINESTREBATTITO_ECG(count,1:length(finestrabattito_ECG)) = finestrabattito_ECG;
%         FINESTREBATTITO_SCG(count,1:length(finestrabattito_SCG)) = finestrabattito_SCG;
%         
%         length_SCG = length(finestrabattito_SCG);
%         length_ECG = length(finestrabattito_ECG);
%         if length_ECG > maxlength_ECG
%             maxlength_ECG = length_ECG
%         end 
%         if length_SCG > maxlength_SCG
%             maxlength_SCG = length_SCG
%         end 
  %%      
        [pks,locs] = findpeaks(finestrabattito_SCG);
        [pksNeg,locsNeg] = findpeaks(-finestrabattito_SCG);
        pksNeg = - pksNeg;

        % Voglio calcolare i primi 3 massimi sapendo che sono dopo il picco
        % R --> Calcolo la posizione all'interno della finestra del picco
        locsqrs1 = round(window_SCG);
        locsafterpicco = find(locs>locsqrs1);
        for p = 1:length(locsafterpicco)
            locsafterR(p) = locs(locsafterpicco(p));
            pksafterR(p) = pks(locsafterpicco(p));
        end 
        %%
        if length(locsafterpicco)>=3
            count = count+1;
            FINESTREBATTITO_ECG(count,1:length(finestrabattito_ECG)) = finestrabattito_ECG;
            FINESTREBATTITO_SCG(count,1:length(finestrabattito_SCG)) = finestrabattito_SCG;
            
            length_SCG = length(finestrabattito_SCG);
            length_ECG = length(finestrabattito_ECG);
            if length_ECG > maxlength_ECG
                maxlength_ECG = length_ECG
            end 
            if length_SCG > maxlength_SCG
                maxlength_SCG = length_SCG
            end 
        
        %%
        %% AO: aortic valve opening
            maxpks = maxk(pksafterR,3) % ho i 3 picchi di ampiezza maggiore dopo R
            for p = 1:length(maxpks)
                maxlocs(p) = find(finestrabattito_SCG == maxpks(p));
                maxp(p,:) = [maxpks(p),maxlocs(p)];
            end
            maxsort = sortrows(maxp,2);
            AO(i,:) = [qrs1-window_SCG+maxsort(1,2)-1 maxsort(1,1)];
        %% RE: rapid ejection
            poslocsRE = find(locsafterR > maxsort(1,2));
            locsRE = locsafterR(poslocsRE(1)); %primo picco dopo AO
            posRE = pksafterR(poslocsRE(1));
            RE(i,:) = [qrs1-window_SCG+locsRE-1 posRE];
        %% T waves and AC: aortic valve closure
     % INDIVIDUO MASSIMO DI ONDA T, il picco AC si trova dopo
           [picchiECG,locspicchiECG] = findpeaks(finestrabattito_ECG); 
           piccoR = [find(finestrabattito_ECG == max(picchiECG)) max(picchiECG)]
           locsdopoR_ECG = find(locspicchiECG > piccoR(1))
           pksdopoR_ECG = picchiECG(locsdopoR_ECG)
           % POSSO CONTARE IL NUMERO DI POSITIVI DOPO ECG
           numero = 0;
           for p = 1:length(pksdopoR_ECG)
               if (pksdopoR_ECG(p) > 0)
                   % posso mettere una condizione sul tempo
                   % tipo che siano lontani più di 0,05s -> 0.05*1024 = 51,2 
                   numero = numero + 1
               end 
           end 
           if numero == 1 %prendo il massimo, sono certa di prendere solo picchi positivi e non un picco positivo ed il minor negativo
            pksT = max(pksdopoR_ECG);
            locs_T = find(finestrabattito_ECG == pksT);
           else % se ho 2 picchi positivi, prendo il primo in ordine temporale 
            pks2T = maxk(pksdopoR_ECG,2);
            for p = 1:length(pks2T)
                locs_2T(p) = find(finestrabattito_ECG == pks2T(p));
                dueT(p,:) = [pks2T(p) locs_2T(p)];
            end
            dueT = sortrows(dueT,2);
            locs_T = dueT(1,2);
            pksT = dueT(1,1);
           end 
           % prendo T come il primo picco in ordine temporale, non è detto
           % infatti che sia per forza il picco maggiore
           T(i,:) = [qrs_I(i)-window_ECG+locs_T-1 pksT];
    %       afterT = finestrabattito_ECG(locs_T:length(finestrabattito_ECG))
          % cerco il minimo dopo l'onda T -> punto in cui più o meno torna ad
          % isoelettrica
           [pksNeg_ECG,locsNeg_ECG] = findpeaks(-finestrabattito_ECG);
           pksNeg_ECG = - pksNeg_ECG;
           locsmindopoT = find(locsNeg_ECG > locs_T);
           pksmindopoT = pksNeg_ECG(locsmindopoT);
           % vogliamo che i picchi siano solo negativi, di quelli negativi
           % prendiamo poi il primo
           negativi = find(pksmindopoT<0);
           if ~isempty(negativi)
               pksmindopoT = pksmindopoT(negativi(1):negativi(end));
               locsmindopoT = locsmindopoT(negativi(1):negativi(end));
           end 
           %% se locsmindopoT è vuoto, non calcolo la fine dell'onda T e prendo direttamente il primo massimo dopo RE
           if isempty(locsmindopoT) % 
               posAC = max(pksafterR(poslocsRE(2:end)));
               posizioneAC = find(finestrabattito_SCG == posAC);
               AC(i,:) = [qrs1-window_SCG+posizioneAC-1 posAC];
           else % se posso calcolare fine onda T
               fine_T = [qrs_I(i)-window_ECG+locsNeg_ECG(locsmindopoT(1))-1 pksmindopoT(1)];
               fine_T_conv = (fine_T(1)/1024)*64;
               fine_T_SCG = fine_T_conv-qrs1+window_SCG+1;
               locsmaxdopofineT = find(locs > fine_T_SCG);
               pksmaxdopofineT = pks(locsmaxdopofineT);
               possibileAC = maxk(pksmaxdopofineT,2);
               if isempty(possibileAC) || (possibileAC(1)<0)
                   % prendo il picco appena prima della fine dell'onda T
                   locsprimadifineT = find(locs < fine_T_SCG);
                   nuovalocsAC = locs(locsprimadifineT(end));
                   nuovapksAC = pks(locsprimadifineT(end));
                   AC(i,:) = [qrs1-window_SCG+nuovalocsAC-1 nuovapksAC];
               else
                   AC(i,:) = [qrs1-window_SCG+find(finestrabattito_SCG == possibileAC(1))-1 possibileAC(1)];
               end
           end 

       %% MC: mitral valve closure 
       % è il picco prima di AO
            locsbeforeAO = find(locs < maxsort(1,2));
            MClocs = locs(locsbeforeAO(end));
            MCpks = pks(locsbeforeAO(end));
            MC(i,:) = [qrs1-window_SCG+MClocs-1 MCpks];
       %% IVC: isovolumteric contraction 
            for p = 1:length(pksNeg)
                if (locsNeg(p) < maxsort(1,2) && locsNeg(p) > MClocs)
                    peakIVC = [pksNeg(p),locsNeg(p)];
                end
            end 
            if peakIVC(2) > locsqrs1
                IVC(i,:) = [qrs1-window_SCG+peakIVC(2)-1 peakIVC(1)];
            else 
                % IVC primo minimo dopo R, AO e' poi il picco successivo , RE il secondo ed MC è
                % il picco prima. Non cambia niente per AC.
                locs_negdopoR = find(locsNeg > locsqrs1);
                pks_negdopoR = pksNeg(locs_negdopoR(1));
                new_posIVC = locsNeg(locs_negdopoR(1));
                IVC(i,:) = [qrs1-window_SCG+new_posIVC-1 pks_negdopoR];
    
                locsafterIVC = find(locs > new_posIVC);
                AOlocs = locs(locsafterIVC(1));
                AOpks = pks(locsafterIVC(1));
                AO(i,:) = [qrs1-window_SCG+AOlocs-1 AOpks];
                RE(i,:) = [qrs1-window_SCG+locs(locsafterIVC(2))-1 pks(locsafterIVC(2))];
                locsbeforeIVC = find(locs < new_posIVC);
                MC(i,:) = [qrs1-window_SCG+locs(locsbeforeIVC(end))-1 pks(locsbeforeIVC(end))];
                % se AC corrisponde con RE, dico che AC e' il picco
                % successivo... ha senso? direi di si ma non sono certa di
                % questa cosa. A noi non interessa l'onda RE ma c'è comunque un
                % ordine temporale ed è presente!
                if (RE(i,1) == AC(i,1))
                    if (length(locsafterIVC)>=3)
                        AC(i,:) = [qrs1-window_SCG+locs(locsafterIVC(3))-1 pks(locsafterIVC(3))];
                    end 
                end 
                % e se non ci fosse un terzo picco...? allora rimane così
                % com'è... FORSE 
            end
        %% check su RE ed AC
            if (RE(i,1) == AC(i,1)) 
                if (length(poslocsRE)>=2)
                    locsACnew = locsafterR(poslocsRE(2)); %primo picco dopo AO
                    posACnew = pksafterR(poslocsRE(2));
                    AC(i,:) = [qrs1-window_SCG+locsACnew-1 posACnew];
                end 
            end 
%% 
            picchi_parziali(count,1) = n_picchi;
        %%
%         figure()
%         a = subplot(211)
%         plot((qrs_I(i)-window_ECG:qrs_I(i+1)-window_ECG)./1024,finestrabattito_ECG),xlabel('[s]'); hold on; plot(qrs_I(i)./1024,qrs_AMP(i),'*r'); hold on; plot(T(i,1)/1024,T(i,2),'*g'); hold on;
%         xline(qrs_I(i)/1024); hold on; xline(T(i,1)/1024,'--g'); hold on;
%         %plot(fine_T(1)/1024,fine_T(2),'*b'); hold on; %xline(fine_T(1)/1024,'--b')
%         b = subplot(212)
%         plot((qrs1-window_SCG:qrs2-window_SCG)./64,finestrabattito_SCG),xlabel('[s]'); hold on; 
%         for r = 1:length(row)
%             plot((POS_picchi_SCG(row(r)))./64,AMP_picchi_SCG(row(r)),'mo')
%         end 
%         xline(qrs1/64); hold on; xline(T(i,1)/1024,'--g'); hold on; %xline(fine_T(1)/1024,'--b')
%         xline((AO(i,1)+window_SCG)./64); hold on
%         plot(IVC(i,1)/64,IVC(i,2),'*r'); hold on
%         plot(AO(i,1)/64,AO(i,2),'*r'); hold on;
%         plot(RE(i,1)/64,RE(i,2),'*r'); hold on;
%         plot(AC(i,1)/64,AC(i,2),'*r'); hold on;
%         plot(MC(i,1)/64,MC(i,2),'*r'); hold on;
%         line([AO(i,1)/64 IVC(i,1)/64],[AO(i,2) IVC(i,2)])
%         text(IVC(i,1)/64,IVC(i,2),' IVC')
%         text(AO(i,1)/64,AO(i,2),' AO')
%         text(RE(i,1)/64,RE(i,2),' RE')
%         text(AC(i,1)/64,AC(i,2),' AC')
%         text(MC(i,1)/64,MC(i,2),' MC')
%         sgtitle(i)
%          pause

%% Estraggo i parameters 
            [tIVCAO,tIVCAC,ampIVCAO,ampIVCAC,slopeIVCAO] = extractfeatures(AO(i,1),AO(i,2),IVC(i,1),IVC(i,2),AC(i,1),AC(i,2),64);
            t_IVCAO(i,1) = tIVCAO;
            t_IVCAC(i,1) = tIVCAC;
            amp_IVCAO(i,1) = ampIVCAO;
            amp_IVCAC(i,1) = ampIVCAC;
            if isnan(slopeIVCAO)
                slope_IVCAO(i,1) = 0;
            else
                slope_IVCAO(i,1) = slopeIVCAO;
            end   
        end % chiude il check dei 3 picchi dopo R
    end % chiude il numero di picchi = 1 o 2

    picchi_totali(i,1) = n_picchi;
    

    clearvars locsafterpicco locsafterR pksafterR maxpks maxlocs maxp maxsort ...
        poslocsRE locsdopoR_ECG pksdopoR_ECG numero locsmindopoT pksmindopoT ...
        locsmaxdopofineT pksmaxdopofineT locsprimadifineT possibileAC locsbeforeAO ...
        peakIVC locs_negdopoR pks_negdopoR locsafterIVC locsbeforeIVC pksmindopoT ...
        locsmindopoT
 end

    %Tolgo le colonne che non avevo messo a zero a caso
    FINESTREBATTITO_SCG = FINESTREBATTITO_SCG(:,1:maxlength_SCG);
    FINESTREBATTITO_ECG = FINESTREBATTITO_ECG(:,1:maxlength_ECG);

    % Prima di salvare tolgo tutte le righe prima di iniziopicchi!!!
    % IMPORTANTE!!!! [] devo eliminarle da tutti
    MC(1:iniziopicchi-1,:) = []; 
    RE(1:iniziopicchi-1,:) = []; 
    AO(1:iniziopicchi-1,:) = []; 
    AC(1:iniziopicchi-1,:) = []; 
    IVC(1:iniziopicchi-1,:) = []; 
    t_IVCAO(1:iniziopicchi-1,:) = []; 
    t_IVCAC(1:iniziopicchi-1,:) = [];
    amp_IVCAO(1:iniziopicchi-1,:) = [];
    amp_IVCAC(1:iniziopicchi-1,:) = [];
    slope_IVCAO(1:iniziopicchi-1,:) = [];
    picchi_totali(1:iniziopicchi-1,:) = [];
    FINESTREBATTITO_ECG(count+1:end,:) = [];
    FINESTREBATTITO_SCG(count+1:end,:) = [];
    R(1:iniziopicchi-1,:) = [];cd 

    % Salvo i dati (fiducial points e parameters)
    name = erase(name,"ECG_FILT-")
    save(['C:\Users\feder\Desktop\Tesi\Data\Parameters SCG\' 'Parameters SCG-' name],'t_IVCAO','t_IVCAC','amp_IVCAO','amp_IVCAC','slope_IVCAO')
    save(['C:\Users\feder\Desktop\Tesi\Data\Fiducial Points SCG\' 'Fiducials SCG-' name],'AO','RE','IVC','AC','MC','R')
    save(['C:\Users\feder\Desktop\Tesi\Data\Windows\' 'Windows-' name],'picchi_totali','picchi_parziali','FINESTREBATTITO_ECG','FINESTREBATTITO_SCG')

%% Variable 'FINESTREBATTITO_ECG' was not saved. For variables larger than 2GB use MAT-file version 7.3 or later --> PROBLEMA NEL SALVATAGGIO DI FINESTRANATTITO_ECG  
end % chiude il numero di soggetti



