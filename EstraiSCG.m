%% Define fiducial points and extract paramerers - SCG 
%% Qui voglio estrarre l'SCG solo se è presente anche l'ECG.
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
 
% Fiducial Points
    MC = zeros(peaksFORwindow(end-1),2); % in realtà non mi serve
    RF = zeros(peaksFORwindow(end-1),2); % in realtà non mi serve
    AO = zeros(peaksFORwindow(end-1),2);
    AC = zeros(peaksFORwindow(end-1),2);
    IVC = zeros(peaksFORwindow(end-1),2); 
% Parameters from SCG
    t_IVCAO = zeros(peaksFORwindow(end-1),2);
    t_IVCAC = zeros(peaksFORwindow(end-1),2);
    amp_IVCAO = zeros(peaksFORwindow(end-1),2);
    amp_IVCAC = zeros(peaksFORwindow(end-1),2);
    slope_IVCAO = zeros(peaksFORwindow(end-1),2);

    %% non prendo dal primo picco - taglio i primi 5 minuti 
    cinqueminuti = 150;
    iniziopicchi = peaksFORwindow(cinqueminuti);

 %%
 for i = iniziopicchi:peaksFORwindow(end-1) %scorro tutti i picchi ECG
% for i = 400:420
    qrs1 = (qrs_I(i)/1024)*64;
    qrs2 = (qrs_I(i+1)/1024)*64;
    finestra_dopo = qrs2
    finestra_prima = qrs1-window_SCG
    [row]=find(POS_picchi_SCG<finestra_dopo & POS_picchi_SCG>finestra_prima)
    n_picchi = size(row,1)
    
    if (n_picchi == 0)
        % non faccio niente
    elseif (n_picchi >= 3)
        % non faccio niente, probabilmente il segnale è troppo corrotto 
    else 
%     elseif (n_picchi == 1)
        picchi = POS_picchi_SCG(row)
        finestrabattito_ECG = ECG_filt(qrs_I(i)-window_ECG:qrs_I(i+1))';
        finestrabattito_SCG = Acc_z(qrs1-window_SCG:qrs2)';

        if i == iniziopicchi 
        FINESTREBATTITO_ECG = {finestrabattito_ECG};
        FINESTREBATTITO_SCG = {finestrabattito_SCG};
        else 
            FINESTREBATTITO_ECG = {FINESTREBATTITO_ECG,finestrabattito_ECG};
            FINESTREBATTITO_SCG = {FINESTREBATTITO_SCG,finestrabattito_SCG};
        end 
        % All_variables = {a,b,c}
  %%      
        [pks,locs] = findpeaks(finestrabattito_SCG)
        [pksNeg,locsNeg] = findpeaks(-finestrabattito_SCG)
        pksNeg = - pksNeg;

        % Voglio calcolare i primi 3 massimi sapendo che sono dopo il picco
        % R --> Calcolo la posizione all'interno della finestra del picco
        locsqrs1 = round(window_SCG);% sarebbe 12,8 ma le posizioni sono sempre interi, meglio mettere round 

        locsafterpicco = find(locs>locsqrs1)
        for p = 1:length(locsafterpicco)
            locsafterR(p) = locs(locsafterpicco(p))
            pksafterR(p) = pks(locsafterpicco(p))
        end 

        maxpks = maxk(pksafterR,3) % ho i 3 picchi di ampiezza maggiore dopo R
        for p = 1:length(maxpks)
            maxlocs(p) = find(finestrabattito_SCG == maxpks(p))
            maxp(p,:) = [maxpks(p),maxlocs(p)]
        end
        maxsort = sortrows(maxp,2);
        AO(i,:) = [qrs1-window_SCG+maxsort(1,2)-1 maxsort(1,1)];
        RF(i,:) = [qrs1-window_SCG+maxsort(2,2)-1 maxsort(2,1)];

        possibleXAC = qrs1-window_SCG+maxsort(3,2)-1;
        distance = (possibleXAC - AO(i,1))./64;
        if distance >= 0.2 %distanza maggiore di 200 ms, è il terzo massimo 
            AC(i,:) = [qrs1-window_SCG+maxsort(3,2)-1 maxsort(3,1)];
        else % vado dopo AO + 200 ms e cerco il massimo
            locsafter200ms = find(locs>(maxsort(1,2)+window_SCG)) %locs>AO + 200 ms
            for p = 1:length(locsafter200ms)
                locspossibleAC(p) = locs(locsafter200ms(p))
                pkspossibleAC(p) = pks(locsafter200ms(p))
            end 
%             pksAC = max(pkspossibleAC);
%             locsAC = find(finestrabattito_SCG == pksAC);
%             AC(i,:) = [qrs1-window_SCG+locsAC-1 pksAC];
            % PRENDO I 3 MASSIMI E CONSIDERO IL PRIMO IN ORDINE TEMPORALE 
            pksAC = maxk(pkspossibleAC,3);
            for p = 1:length(pksAC)
                locsAC(p) = find(finestrabattito_SCG == pksAC(p))
                pkslocsAC(p,:) = [pksAC(p),locsAC(p)]
            end
            sortAC = sortrows(pkslocsAC)
            AC(i,:) = [qrs1-window_SCG+sortAC(1,2)-1 sortAC(1,1)];
%             AC(i,:) = [qrs1-window_SCG+locsAC-1 pksAC]
        end 


        % MC, è il picco prima di AO
        locsbeforeAO = find(locs < maxsort(1,2))
        MClocs = locs(locsbeforeAO(end))
        MCpks = pks(locsbeforeAO(end))
        MC(i,:) = [qrs1-window_SCG+MClocs-1 MCpks];

        
     
%         minpks = mink(pksNeg,3) % ho i 3 picchi negativi di ampiezza maggiore
%         for p = 1:3
%             minlocs(p) = find(finestrabattito_SCG(i,:) == minpks(p))
%             min(p,:) = [minpks(p),minlocs(p)]
%         end
%         minsort = sortrows(min,2)
%         for p = 1:3
%             if minsort(p,2) < maxsort(1,2)
%                 peakIVC = minsort(p,:)
%             end
%         end

        for p = 1:length(pksNeg)
            if (locsNeg(p) < maxsort(1,2) && locsNeg(p) > MClocs)
                peakIVC = [pksNeg(p),locsNeg(p)]
            end
        end 
       IVC(i,:) = [qrs1-window_SCG+peakIVC(2)-1 peakIVC(1)];
       
        %%
        figure()
        a = subplot(211)
        plot((qrs_I(i)-window_ECG:qrs_I(i+1))./1024,finestrabattito_ECG),xlabel('[s]'); hold on; plot(qrs_I(i)./1024,qrs_AMP(i),'*r'); hold on;
        xline(qrs_I(i)/1024)
        b = subplot(212)
        plot((qrs1-window_SCG:qrs2)./64,finestrabattito_SCG),xlabel('[s]'); hold on; 
        for r = 1:length(row)
            plot((POS_picchi_SCG(row(r)))./64,AMP_picchi_SCG(row(r)),'mo')
        end 
%         plot((POS_picchi_SCG(row(1)):POS_picchi_SCG(row(end)))./64,AMP_picchi_SCG(row(1)):AMP_picchi_SCG(row(end)),'mo'); hold on
        xline(qrs1/64); hold on
        xline((AO(i,1)+window_SCG)./64); hold on
        plot(IVC(i,1)/64,IVC(i,2),'*r'); hold on
        plot(AO(i,1)/64,AO(i,2),'*r'); hold on;
        plot(RF(i,1)/64,RF(i,2),'*r'); hold on;
        %plot(AC(i,1)/64,AC(i,2),'*r'); hold on;
        plot(MC(i,1)/64,MC(i,2),'*r'); hold on;
        line([AO(i,1)/64 IVC(i,1)/64],[AO(i,2) IVC(i,2)])
        text(IVC(i,1)/64,IVC(i,2),' IVC')
        text(AO(i,1)/64,AO(i,2),' AO')
        text(RF(i,1)/64,RF(i,2),' RF')
        %text(AC(i,1)/64,AC(i,2),' AC')
        text(MC(i,1)/64,MC(i,2),' MC')
        sgtitle(i)
%         
%        
%         pause
%%
%     elseif (n_picchi==2)     
%         finestrabattito_ECG = ECG_filt(qrs_I(i)-window_ECG:qrs_I(i+1))';
%         finestrabattito_SCG = Acc_z(qrs1-window_SCG:qrs2)';
%         figure()
%         a = subplot(211)
%         plot((qrs_I(i)-window_ECG:qrs_I(i+1))./1024,finestrabattito_ECG); hold on; plot(qrs_I(i)./1024,qrs_AMP(i),'*r') 
%         b = subplot(212)
%         plot((qrs1-window_SCG:qrs2)./64,finestrabattito_SCG); hold on;
%         for r = 1:length(row)
%             plot((POS_picchi_SCG(row(r)))./64,AMP_picchi_SCG(row(r)),'mo')
%         end 
%         sgtitle(i)
        
        % SE IL PRIMO PICCO E' PRIMA DI R, ALLORA CONSIDERO SOLO IL
        % SECONDO! (i = 200 )
        
%     elseif (n_picchi > 2)
%     else % n_picchi = 0
    end 
%  end 

   
%
%     end
    picchi(i,1) = n_picchi;

    %% QUI METTO LA FUNZIONE CHE ESTRAE I DATI 
    
    [tIVCAO,tIVCAC,ampIVCAO,ampIVCAC,slopeIVCAO] = extractfeatures(AO(i,1),AO(i,2),IVC(i,1),IVC(i,2),AC(i,1),AC(i,2),64)
    t_IVCAO(i,1) = tIVCAO;
    t_IVCAC(i,1) = tIVCAC;
    amp_IVCAO(i,1) = ampIVCAO;
    amp_IVCAC(i,1) = ampIVCAC;
    if isnan(slopeIVCAO)
        slope_IVCAO(i,1) = 0;
    else
        slope_IVCAO(i,1) = slopeIVCAO;
    end 

   
    clearvars locsafterpicco locafterR pksafterR maxpks maxlocs max maxsort locsbeforeAO peakIVC locsafter200ms locspossibileAC pkspossibleAC pksAC locsAC pkslocsAC 

end
    % Salvo i dati (fiducial points e parameters)
    name = erase(name,"ECG-FILT-")
    save(['C:\Users\feder\Desktop\Tesi\Data\Parameters SCG' 'Parameters SCG-' name],'t_IVCAO','t_IVCAC','amp_IVCAO','amp_IVCAC','slope_IVCAO')
    save(['C:\Users\feder\Desktop\Tesi\Data\Fiducial Points SCG' 'Fiducials SCG-' name],'AO','RF','IVC','AC','MC')
%     save(['C:\Users\feder\Desktop\Tesi\Data\Windows SCG' 'Windows SCG-' name],'picchi,FINESTREEEEEEE)

end 



%%
figure()
a = subplot(311)
plot((qrs_I(i)-window_ECG:qrs_I(i+1))./1024,finestrabattito_ECG),xlabel('[s]'); hold on; plot(qrs_I(i)./1024,qrs_AMP(i),'*r'); hold on;
xline(qrs_I(i)/1024)
b = subplot(312)
plot((qrs1-window_SCG:qrs2)./64,finestrabattito_SCG),xlabel('[s]'); hold on; 
for r = 1:length(row)
    plot((POS_picchi_SCG(row(r)))./64,AMP_picchi_SCG(row(r)),'mo')
end 
%         plot((POS_picchi_SCG(row(1)):POS_picchi_SCG(row(end)))./64,AMP_picchi_SCG(row(1)):AMP_picchi_SCG(row(end)),'mo'); hold on
xline(qrs1/64); hold on
%plot(IVC(i,1)/64,IVC(i,2),'*r'); hold on
plot(AO(i,1)/64,AO(i,2),'*r'); hold on;
plot(RF(i,1)/64,RF(i,2),'*r'); hold on;
%plot(AC(i,1)/64,AC(i,2),'*r'); hold on;
%plot(MC(i,1)/64,MC(i,2),'*r'); hold on;
%line([AO(i,1)/64 IVC(i,1)/64],[AO(i,2) IVC(i,2)])+
%text(IVC(i,1)/64,IVC(i,2),' IVC')
text(AO(i,1)/64,AO(i,2),' AO')
text(RF(i,1)/64,RF(i,2),' RF')
%text(AC(i,1)/64,AC(i,2),' AC')
%text(MC(i,1)/64,MC(i,2),' MC')
sgtitle(i)
subplot(313)
plot(finestrabattito_SCG); hold on;xline(locsqrs1)

     