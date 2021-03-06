%% Define fiducial points and extract paramerers - SCG 
%% Qui voglio estrarre l'SCG solo se è presente anche l'ECG.
 clear all
 close all
 clc
%%
% folderPT = 'C:\Users\feder\Desktop\Tesi\Data\PostProc PT 1_10SEC'; 
folderSCG = 'C:\Users\feder\Desktop\Tesi\Data\Picchi SCG - Acc z';
folderPT = 'C:\Users\feder\Desktop\Tesi\Data\PostProc PT 1'; 
folderECG = 'C:\Users\feder\Desktop\Tesi\Data\Filtered ECG';
% folderSCG = 'C:\Users\feder\Desktop\Tesi\Data\Picchi SCG - Acc z - 10 sec';

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
list(end) = [];
N = length(list);
% list(N-1) = [];
% list(N-1) = [];
% N = length(list)

%  addpath 'C:\Users\feder\Desktop\Tesi'\Data\'PostProc PT 1_10SEC'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Filtered ECG'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Data\'PostProc PT 1'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Picchi SCG - Acc z'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Codes\
%  addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Picchi SCG - Acc z - 10 sec'\

% data = xlsread('C:\Users\feder\Desktop\Tesi\Info Pazienti.xlsx','P:P');
[~, sex] = xlsread('C:\Users\feder\Desktop\Tesi\Info Pazienti.xlsx','D:D');
sex(1) = [];
[~,txtdata] = xlsread('C:\Users\feder\Desktop\Tesi\Info Pazienti.xlsx','H:J');
txtdata(1,:) = [];
Inizio_Holter = txtdata(:,1);
Inizio_Periodo_Sonno = txtdata(:,2);
Fine_Periodo_Sonno = txtdata(:,3);

%% Parto da ECG, per ogni battito di ECG considero una finestra a sx e a dx del picco R.. 
for m = 12:12
    
    

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

    % Aggiunto dopo il controllo di RR 
    qrs_I = QRS(:,1)';
    qrs_AMP = QRS(:,2);

    fs_ECG = 1024;
    fs_SCG = 64;
    window_SCG = 0.2*fs_SCG;
    window_ECG = 0.2*fs_ECG;

      iniziopicchi = sum(HR_min(1:5));
      finepicchi = sum(HR_min(end-5:end)); 
      picchi_totali = length(qrs_I)-iniziopicchi-finepicchi; 
      

% Fiducial Points
    MC = zeros(length(qrs_I)-finepicchi,3); % in realtà non mi serve
    RE = zeros(length(qrs_I)-finepicchi,3); % in realtà non mi serve
    AO = zeros(length(qrs_I)-finepicchi,3);
    AC = zeros(length(qrs_I)-finepicchi,3);
    IVC = zeros(length(qrs_I)-finepicchi,3); 
    R = zeros(length(qrs_I)-finepicchi,3);
    Q = zeros(length(qrs_I)-finepicchi,3);
    minAO_RE = zeros(length(qrs_I)-finepicchi,3);
    minbeforeAC = zeros(length(qrs_I)-finepicchi,3);
    T = zeros(length(qrs_I)-finepicchi,3);
    fine_T = zeros(length(qrs_I)-finepicchi,3);
    RR = zeros(length(qrs_I)-finepicchi-1,2);

% Parameters from SCG
    t_IVCAO = zeros(length(qrs_I)-finepicchi,2); t_IVCAC = zeros(length(qrs_I)-finepicchi,2);
    t_IVCMC = zeros(length(qrs_I)-finepicchi,2); t_IVCRE = zeros(length(qrs_I)-finepicchi,2);
    t_IVCminAORE = zeros(length(qrs_I)-finepicchi,2); t_IVCminAC = zeros(length(qrs_I)-finepicchi,2);
    amp_IVCAO = zeros(length(qrs_I)-finepicchi,2); amp_IVCAC = zeros(length(qrs_I)-finepicchi,2);
    amp_IVCMC = zeros(length(qrs_I)-finepicchi,2); amp_IVCRE = zeros(length(qrs_I)-finepicchi,2);
    amp_IVCminAORE = zeros(length(qrs_I)-finepicchi,2); amp_IVCminAC = zeros(length(qrs_I)-finepicchi,2);
    slope_IVCAO = zeros(length(qrs_I)-finepicchi,2); slope_minAORERE = zeros(length(qrs_I)-finepicchi,2);
    slope_minACAC = zeros(length(qrs_I)-finepicchi,2); 
    LVET = zeros(length(qrs_I)-finepicchi,2); % AC-AO
    QS2 = zeros(length(qrs_I)-finepicchi,2); % Q-AC
    QT = zeros(length(qrs_I)-finepicchi,2); % Q-fine onda T
    QTc = zeros(length(qrs_I)-finepicchi,2); % Q-fine onda T (corretto) 
    R_div_T = zeros(length(qrs_I)-finepicchi,2);
    R_AC1 = zeros(length(qrs_I)-finepicchi,2);
    R_MC1 = zeros(length(qrs_I)-finepicchi,2);
    Q_AC1 = zeros(length(qrs_I)-finepicchi,2);
    R_AO = zeros(length(qrs_I)-finepicchi,2);
    
%     FINESTREBATTITO_ECG = zeros(length(qrs_I)-finepicchi,2000);
%     FINESTREBATTITO_SCG = zeros(length(qrs_I)-finepicchi,100);

    maxlength_SCG = 0;
    maxlength_ECG = 0;
    count = 0;
    tag0 = 0; % Analizzo
    tag1 = 0; % 0 picchi 
    tag2 = 0; % < 3 picchi in SCG
    tag3 = 0; % onda T troppo avanti
    tag4 = 0; % no sistole e no diastole
    tag5 = 0; % si sistole, no diastole
    tag6 = 0; % non è presente Q
    tag7 = 0; % finestra troppo corta rispetto a QTc 440 o 460 ms 
    tag8 = 0; % non è presente MC

     % Capisco se siamo durante il giorno o durante la notte
    Inizio = datevec(Inizio_Holter{m});
    Inizio_Sonno = datevec(Inizio_Periodo_Sonno{m});
    Fine_Sonno = datevec(Fine_Periodo_Sonno{m});
    tempofinoasonno_sec = etime(Inizio_Sonno,Inizio); 
    duratasonnodainizio_sec = etime(Fine_Sonno,Inizio);
    tempofinoasonno = tempofinoasonno_sec*fs_ECG;
    duratasonnodainizio = duratasonnodainizio_sec*fs_ECG;
    GN = zeros(length(qrs_I)-finepicchi-1,1);
    n_picchi = zeros(length(qrs_I)-finepicchi-1,1);
 %%
for i = iniziopicchi:length(qrs_I)-finepicchi-1 %scorro tutti i picchi ECG
    if (qrs_I(1,i)>tempofinoasonno && qrs_I(1,i)<duratasonnodainizio)
        GN(i) = 1; % notte
    else
        GN(i) = 0; % giorno
    end
 
    qrs1 = (qrs_I(i)/1024)*64;
    qrs2 = (qrs_I(i+1)/1024)*64;
    finestra_dopo = qrs2-window_SCG; %cambiato
    finestra_prima = qrs1-window_SCG;
    [row]=find(POS_picchi_SCG<finestra_dopo & POS_picchi_SCG>finestra_prima);
    n_picchi(i,1) = size(row,1);
    R(i,1:2) = [qrs_I(i) qrs_AMP(i)];
    RR(i,1:2) = [qrs_I(i+1)-qrs_I(i) qrs_I(i)];

    if (n_picchi(i,1) == 0)
%         tutti i parameters per questa finestra sono NaN
% Si riesce a salvare qualcosa da questa finestra?? 
        tag1 = tag1+1;
        R(i,3) =  1;
        t_IVCAO(i,1) = NaN; t_IVCAC(i,1) = NaN; t_IVCMC(i,1) = NaN; t_IVCRE(i,1) = NaN; t_IVCminAORE(i,1) = NaN; t_IVCminAC(i,1) = NaN;
        amp_IVCAO(i,1) = NaN; amp_IVCAC(i,1) = NaN; amp_IVCMC(i,1) = NaN; amp_IVCRE(i,1) = NaN; amp_IVCminAORE(i,1) = NaN; amp_IVCminAC(i,1) = NaN; 
        slope_IVCAO(i,1) = NaN; slope_minAORERE(i,1) = NaN; slope_minACAC(i,1) = NaN; LVET(i,1) = NaN; QS2(i,1) = NaN; QT(i,1) = NaN; QTc(i,1) = NaN;
        MC(i,1:2) = NaN; RE(i,1:2) = NaN; AO(i,1:2) = NaN; AC(i,1:2) = NaN; IVC(i,1:2) = NaN; minAO_RE(i,1:2) = NaN; minbeforeAC(i,1:2) = NaN;
        Q(i,1:2) = NaN; fine_T(i,1:2) = NaN; T(i,1:2) = NaN; RR(i,1:2) = NaN; R_MC1(i,1) = NaN; R_AC1(i,1) = NaN; R_AO(i,1) = NaN; 
       
    else %% n_picchi == 1 or n_picchi == 2
        picchi = POS_picchi_SCG(row);
        finestrabattito_ECG = ECG_filt(qrs_I(i)-window_ECG:qrs_I(i+1)-window_ECG)';
        finestrabattito_SCG = Acc_z(qrs1-window_SCG:qrs2-window_SCG)';

        [pks,locs] = findpeaks(finestrabattito_SCG);
        [pksNeg,locsNeg] = findpeaks(-finestrabattito_SCG);
        pksNeg = - pksNeg;
        % Voglio calcolare i primi 3 massimi sapendo che sono dopo il picco
        % R --> Calcolo la posizione all'interno della finestra del picco
        locsqrs1 = round(window_SCG);
        locsafterpicco = find(locs>locsqrs1);
        locsafterR = locs(locsafterpicco);
        pksafterR = pks(locsafterpicco);
        % Voglio che nella finestra siano presenti almeno 3 picchi:
        % AO,RE,AC
        if length(locsafterpicco)>=3 
%             count = count+1;
%             FINESTREBATTITO_ECG(count,1:length(finestrabattito_ECG)) = finestrabattito_ECG;
%             FINESTREBATTITO_SCG(count,1:length(finestrabattito_SCG)) = finestrabattito_SCG;
%             
%             length_SCG = length(finestrabattito_SCG);
%             length_ECG = length(finestrabattito_ECG);
%             if length_ECG > maxlength_ECG
%                 maxlength_ECG = length_ECG;
%             end 
%             if length_SCG > maxlength_SCG
%                 maxlength_SCG = length_SCG;
%             end 
        
            %% Cerco il punto Q e poi di conseguenza 
            [picchiECG,locspicchiECG] = findpeaks(finestrabattito_ECG); 
            piccoR = [find(finestrabattito_ECG == max(picchiECG)) max(picchiECG)];
            locsprimaR_ECG = find(locspicchiECG < piccoR(1));
            pksprimaR_ECG = picchiECG(locsprimaR_ECG);
            if isempty(locsprimaR_ECG)
                tag6 = tag6+1;
                R(i,3) = 6;
                t_IVCAO(i,1) = NaN; t_IVCAC(i,1) = NaN; t_IVCMC(i,1) = NaN; t_IVCRE(i,1) = NaN; t_IVCminAORE(i,1) = NaN; t_IVCminAC(i,1) = NaN;
                amp_IVCAO(i,1) = NaN; amp_IVCAC(i,1) = NaN; amp_IVCMC(i,1) = NaN; amp_IVCRE(i,1) = NaN; amp_IVCminAORE(i,1) = NaN; amp_IVCminAC(i,1) = NaN; 
                slope_IVCAO(i,1) = NaN; slope_minAORERE(i,1) = NaN; slope_minACAC(i,1) = NaN; LVET(i,1) = NaN; QS2(i,1) = NaN; QT(i,1) = NaN; QTc(i,1) = NaN;
                MC(i,1:2) = NaN; RE(i,1:2) = NaN; AO(i,1:2) = NaN; AC(i,1:2) = NaN; IVC(i,1:2) = NaN; minAO_RE(i,1:2) = NaN; minbeforeAC(i,1:2) = NaN;
                Q(i,1:2) = NaN; fine_T(i,1:2) = NaN; T(i,1:2) = NaN; RR(i,1:2) = NaN; R_MC1(i,1) = NaN; R_AC1(i,1) = NaN; R_div_T(i,1) = NaN; R_AO(i,1) = NaN;
               continue
            end 

            Q(i,1:2) = [qrs_I(i)-window_ECG+locspicchiECG(locsprimaR_ECG(end))-1 pksprimaR_ECG(end)]; %così è pos in sec in ECG
            Q_SCG_locs = (locspicchiECG(locsprimaR_ECG(end))/1024)*64;
            Q_ECG_locs = (locspicchiECG(locsprimaR_ECG(end)));
            Q_SCG = [qrs1-window_SCG+Q_SCG_locs-1]; 

            % Ricerco il valore corretto di HR in base alla finestra di 5
            % min in cui si trova -> mi serve per tag 5
            qrs_I_min = (qrs_I(i)./fs_ECG)./60;
            ricercariga = round(round(qrs_I_min)/5);
            bpm_min = HR_5min(ricercariga+1);

            if sex{m} == 'F'
               QS2maxline = -0.0020*bpm_min+0.549;
           else
               QS2maxline = -0.0021*bpm_min+0.546;
            end 
            QS2maxlined = qrs1-window_SCG+(Q_SCG_locs+QS2maxline*64)-1;
            %% AO: aortic valve opening
            maxpks = maxk(pksafterR,3); % ho i 3 picchi di ampiezza maggiore dopo R
            p = 0;
            for p = 1:length(maxpks)
                maxlocs(p) = find(finestrabattito_SCG == maxpks(p));
                maxp(p,:) = [maxpks(p),maxlocs(p)];
            end
            maxsort = sortrows(maxp,2);
            AO(i,1:2) = [qrs1-window_SCG+maxsort(1,2)-1 maxsort(1,1)];
            %% RE: rapid ejection
            poslocsRE = find(locsafterR > maxsort(1,2));
            locsRE = locsafterR(poslocsRE(1)); %primo picco dopo AO
            posRE = pksafterR(poslocsRE(1));
            RE(i,1:2) = [qrs1-window_SCG+locsRE-1 posRE];
            %% T waves and AC: aortic valve closure
         % INDIVIDUO MASSIMO DI ONDA T, il picco AC si trova dopo

           locsdopoR_ECG = locspicchiECG(find(locspicchiECG > piccoR(1)));
%              pksdopoR_ECG = picchiECG(locsdopoR_ECG)
           pksdopoR_ECG = picchiECG(find(locspicchiECG > piccoR(1)));

           % POSSO CONTARE IL NUMERO DI POSITIVI DOPO ECG
           numero = 0;
           p = 0;
           for p = 1:length(pksdopoR_ECG)
               if (pksdopoR_ECG(p) > 0)
                   numero = numero + 1;
               end 
           end 
           if numero == 1 %prendo il massimo, sono certa di prendere solo picchi positivi e non un picco positivo ed il minor negativo
               pksT = max(pksdopoR_ECG);
               locs_T = find(finestrabattito_ECG == pksT);
           else % se ho 2 picchi positivi, prendo il primo in ordine temporale 
                pks2T = maxk(pksdopoR_ECG,2);
                p = 0;
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

           if sex{m} == 'F'
               QTmax_cost = 460;
           else
               QTmax_cost = 440;
           end 
            RR_sys = (40*RR(i))./100; % valore a partire da R
            RR_sys = RR_sys+window_ECG;
            RR_sysd = qrs_I(i)-window_ECG+RR_sys-1;
            
            if  locs_T> RR_sys
               % La posizione di T è troppo avanti -> elimino tutto
               tag3 = tag3+1;
               R(i,3) = 3;
               t_IVCAO(i,1) = NaN; t_IVCAC(i,1) = NaN; t_IVCMC(i,1) = NaN; t_IVCRE(i,1) = NaN; t_IVCminAORE(i,1) = NaN; t_IVCminAC(i,1) = NaN;
               amp_IVCAO(i,1) = NaN; amp_IVCAC(i,1) = NaN; amp_IVCMC(i,1) = NaN; amp_IVCRE(i,1) = NaN; amp_IVCminAORE(i,1) = NaN; amp_IVCminAC(i,1) = NaN; 
               slope_IVCAO(i,1) = NaN; slope_minAORERE(i,1) = NaN; slope_minACAC(i,1) = NaN; LVET(i,1) = NaN; QS2(i,1) = NaN; QT(i,1) = NaN; QTc(i,1) = NaN;
               MC(i,1:2) = NaN; RE(i,1:2) = NaN; AO(i,1:2) = NaN; AC(i,1:2) = NaN; IVC(i,1:2) = NaN; minAO_RE(i,1:2) = NaN; minbeforeAC(i,1:2) = NaN;
               Q(i,1:2) = NaN; fine_T(i,1:2) = NaN; T(i,1:2) = NaN; RR(i,1:2) = NaN; R_MC1(i,1) = NaN; R_AC1(i,1) = NaN; R_div_T(i,1) = NaN; R_AO(i,1) = NaN;
               continue
            end 
            if (length(finestrabattito_ECG)-Q_ECG_locs)<QTmax_cost
                % la mia finestra finisce prima della costante, è una
                % finestra troppo corta --> forse potrei studiarla in un
                % altro modo ma al momento non lo so
                tag7 = tag7+1;
                R(i,3) = 7;
                t_IVCAO(i,1) = NaN; t_IVCAC(i,1) = NaN; t_IVCMC(i,1) = NaN; t_IVCRE(i,1) = NaN; t_IVCminAORE(i,1) = NaN; t_IVCminAC(i,1) = NaN;
                amp_IVCAO(i,1) = NaN; amp_IVCAC(i,1) = NaN; amp_IVCMC(i,1) = NaN; amp_IVCRE(i,1) = NaN; amp_IVCminAORE(i,1) = NaN; amp_IVCminAC(i,1) = NaN; 
                slope_IVCAO(i,1) = NaN; slope_minAORERE(i,1) = NaN; slope_minACAC(i,1) = NaN; LVET(i,1) = NaN; QS2(i,1) = NaN; QT(i,1) = NaN; QTc(i,1) = NaN;
                MC(i,1:2) = NaN; RE(i,1:2) = NaN; AO(i,1:2) = NaN; AC(i,1:2) = NaN; IVC(i,1:2) = NaN; minAO_RE(i,1:2) = NaN; minbeforeAC(i,1:2) = NaN;
                Q(i,1:2) = NaN; fine_T(i,1:2) = NaN; T(i,1:2) = NaN; RR(i,1:2) = NaN; R_MC1(i,1) = NaN; R_AC1(i,1) = NaN; R_div_T(i,1) = NaN; R_AO(i,1) = NaN;
                continue
            end 

                count = count+1;
%                 FINESTREBATTITO_ECG(count,1:length(finestrabattito_ECG)) = finestrabattito_ECG;
%                 FINESTREBATTITO_SCG(count,1:length(finestrabattito_SCG)) = finestrabattito_SCG;
                length_SCG = length(finestrabattito_SCG);
                length_ECG = length(finestrabattito_ECG);
                if length_ECG > maxlength_ECG
                    maxlength_ECG = length_ECG;
                end 
                if length_SCG > maxlength_SCG
                    maxlength_SCG = length_SCG;
                end 

               T(i,1:2) = [qrs_I(i)-window_ECG+locs_T-1 pksT];
               % Cerco la fine dell'onda T
               T40 = locs_T + 0.04*1024;
               T40d = qrs_I(i)-window_ECG+T40-1;
               T80 = locs_T + 0.08*1024;
               T80d = qrs_I(i)-window_ECG+T80-1;
               % se T è abbastanza avanti ed il segnale è corto, devo ridurre la finestra di xr
               if length(finestrabattito_ECG)<T80
                   T80 = length(finestrabattito_ECG);
                   T40 = (T80-locs_T)/2+locs_T;
               end 
    
               finestraxm = finestrabattito_ECG(locs_T:T40);
               derivatam = gradient(finestraxm);
               xm = find(derivatam == max(derivatam));
               ym = finestraxm(xm);
               Xm = locs_T+xm;
               Xmd = qrs_I(i)-window_ECG+Xm-1;
    
               finestraxr = finestrabattito_ECG(T40:T80);
               derivatar = gradient(finestraxr);
               absderivatar = abs(derivatar);
               xr = find(absderivatar == min(absderivatar));
               yr = finestraxr(xr);
               Xr = T40+xr;
               Xrd = qrs_I(i)-window_ECG+Xr-1;
    
                % mi interessa l'area massima 
               Area = 0;
               Area1 = zeros(round(Xr)-round(Xm),1);
               step = 0;
               for x = round(Xm):round(Xr)
                   step = step+1;
                   y = finestrabattito_ECG(x);
                   A = 1/2*(ym-y)*(2*Xr-Xm-x);
                   Area1(step) = A;
                   if abs(A) > Area
                       Area = abs(A);
                       xi = x;
                   end 
               end 
               yi = finestrabattito_ECG(xi);
               Xid = qrs_I(i)-window_ECG+xi-1;
    
               fine_T(i,1:2) = [Xid yi];
               fine_T_conv = (fine_T(i,1)/1024)*64;
               fine_T_SCG = fine_T_conv-qrs1+window_SCG+1;
               locsmaxdopofineT = find(locs > fine_T_SCG);
               pksmaxdopofineT = pks(locsmaxdopofineT);
               if ~isempty(locsmaxdopofineT) % se esiste allora posso trovare AC ed il minimo prima di AC
                   nuovalocsAC = locs(locsmaxdopofineT(1));
                   nuovapksAC = pks(locsmaxdopofineT(1));
                   AC(i,1:2) = [qrs1-window_SCG+nuovalocsAC-1 nuovapksAC];
                    %% Picco minimo prima di AC (dal punto di vista temporale)
                   LOCSAC = AC(i,1)-qrs1+window_SCG+1;
                   locsbeforeAC = find(locsNeg < LOCSAC);
                   minbeforeAC(i,1:2) = [qrs1-window_SCG+locsNeg(locsbeforeAC(end))-1 pksNeg(locsbeforeAC(end))];
                    %% check su RE ed AC
                    if (RE(i,1) == AC(i,1)) 
                        if (length(poslocsRE)>=2)
                            locsACnew = locsafterR(poslocsRE(2)); %primo picco dopo AO
                            posACnew = pksafterR(poslocsRE(2));
                            AC(i,1:2) = [qrs1-window_SCG+locsACnew-1 posACnew];
                        end 
                    end 
               end 
      
           %% MC: mitral valve closure 
           % è il picco prima di AO
            locsbeforeAO = find(locs < maxsort(1,2));
            if isempty(locsbeforeAO)
                tag8 = tag8+1;
                R(i,3) = 8;
                t_IVCAO(i,1) = NaN; t_IVCAC(i,1) = NaN; t_IVCMC(i,1) = NaN; t_IVCRE(i,1) = NaN; t_IVCminAORE(i,1) = NaN; t_IVCminAC(i,1) = NaN;
                amp_IVCAO(i,1) = NaN; amp_IVCAC(i,1) = NaN; amp_IVCMC(i,1) = NaN; amp_IVCRE(i,1) = NaN; amp_IVCminAORE(i,1) = NaN; amp_IVCminAC(i,1) = NaN; 
                slope_IVCAO(i,1) = NaN; slope_minAORERE(i,1) = NaN; slope_minACAC(i,1) = NaN; LVET(i,1) = NaN; QS2(i,1) = NaN; QT(i,1) = NaN; QTc(i,1) = NaN;
                MC(i,1:2) = NaN; RE(i,1:2) = NaN; AO(i,1:2) = NaN; AC(i,1:2) = NaN; IVC(i,1:2) = NaN; minAO_RE(i,1:2) = NaN; minbeforeAC(i,1:2) = NaN;
                Q(i,1:2) = NaN; fine_T(i,1:2) = NaN; T(i,1:2) = NaN; RR(i,1:2) = NaN; R_MC1(i,1) = NaN; R_AC1(i,1) = NaN; R_div_T(i,1) = NaN; R_AO(i,1) = NaN;
                continue
            end 
            MClocs = locs(locsbeforeAO(end));
            MCpks = pks(locsbeforeAO(end));
            MC(i,1:2) = [qrs1-window_SCG+MClocs-1 MCpks];
           %% IVC: isovolumteric contraction 
           p = 0;
            for p = 1:length(pksNeg)
                if (locsNeg(p) < maxsort(1,2) && locsNeg(p) > MClocs)
                    peakIVC = [pksNeg(p),locsNeg(p)];
                end
            end 
            if peakIVC(2) > locsqrs1
                IVC(i,1:2) = [qrs1-window_SCG+peakIVC(2)-1 peakIVC(1)];
            else 
                % IVC primo minimo dopo R, AO e' poi il picco successivo , RE il secondo ed MC è
                % il picco prima. Non cambia niente per AC.
                locs_negdopoR = find(locsNeg > locsqrs1);
                pks_negdopoR = pksNeg(locs_negdopoR(1));
                new_posIVC = locsNeg(locs_negdopoR(1));
                IVC(i,1:2) = [qrs1-window_SCG+new_posIVC-1 pks_negdopoR];
    
                locsafterIVC = find(locs > new_posIVC);
                AOlocs = locs(locsafterIVC(1));
                AOpks = pks(locsafterIVC(1));
                AO(i,1:2) = [qrs1-window_SCG+AOlocs-1 AOpks];
                RE(i,1:2) = [qrs1-window_SCG+locs(locsafterIVC(2))-1 pks(locsafterIVC(2))];
                locsbeforeIVC = find(locs < new_posIVC);
                MC(i,1:2) = [qrs1-window_SCG+locs(locsbeforeIVC(end))-1 pks(locsbeforeIVC(end))];
                % se AC corrisponde con RE, dico che AC e' il picco
                % successivo... ha senso? direi di si ma non sono certa di
                % questa cosa. A noi non interessa l'onda RE ma c'è comunque un
                % ordine temporale ed è presente!
                if (RE(i,1) == AC(i,1))
                    if (length(locsafterIVC)>=3)
                        AC(i,1:2) = [qrs1-window_SCG+locs(locsafterIVC(3))-1 pks(locsafterIVC(3))];
                    end 
                end 
            end

           %% Picco di minimo tra AO ed RE
           p = 0;
           for p = 1:length(pksNeg)
                if (locsNeg(p) > maxsort(1,2) && locsNeg(p) < locsRE)
                    peakminAO_RE = [pksNeg(p),locsNeg(p)];
                end
           end 
           % se dovesse essercene più di uno, li metto in ordine dal più
           % negativo ad il meno negativo e prendo il più negativo
           if (size(peakminAO_RE,1)>1)
               peakminAO_RE = sortrows(peakminAO_RE,1); % riordino in base al valore della prima colonna, ampiezza dei picchi 
               peakminAO_RE = peakminAO_RE(1,:); % la prima riga dovrebbe corrispondere al valore più piccolo (più negativo)
           end 
           minAO_RE(i,1:2) = [qrs1-window_SCG+peakminAO_RE(2)-1 peakminAO_RE(1)];
         

    %% 
                picchi_parziali(count,1) = n_picchi(i,1);
%                 R_SCG = (R(i,1)/1024)*64;
%                 R_MC1(i,1) = R_SCG-MC(i,1);
%             
            figure()
            a = subplot(211)
            plot((qrs_I(i)-window_ECG:qrs_I(i+1)-window_ECG)./1024,finestrabattito_ECG),xlabel('[s]'); hold on; plot(qrs_I(i)./1024,qrs_AMP(i),'*r'); hold on; plot(T(i,1)/1024,T(i,2),'*g'); hold on;
            xline(qrs_I(i)/1024); hold on; xline(T(i,1)/1024,'--g'); hold on;
            xline(T40d/1024,'--m'); hold on; xline(T80d/1024,'--m'); hold on;
            plot(Xmd/1024,ym,'*m'); hold on; plot(Xrd/1024,yr,'*m'); hold on;
            plot(Xid/1024,yi,'*b'); hold on;
            plot(Q(i,1)/1024,Q(i,2),'*b'); hold on;
%             xline(QS2maxline_ECGd/1024,'--r'); hold on;
              xline(RR_sysd/1024,'--y'); hold on;
            plot(fine_T(i,1)/1024,fine_T(i,2),'*b');
            hold on; xline(fine_T(i,1)/1024,'--b')
            text(qrs_I(i)./1024,qrs_AMP(i),' R')
            text(T(i,1)/1024,T(i,2),' T')
            text(fine_T(1)/1024,fine_T(2),' fine onda T')
            text(Q(i,1)/1024,Q(i,2),' Q')
            b = subplot(212)
            plot((qrs1-window_SCG:qrs2-window_SCG)./64,finestrabattito_SCG),xlabel('[s]'); hold on; 
            for r = 1:length(row)
                plot((POS_picchi_SCG(row(r)))./64,AMP_picchi_SCG(row(r)),'mo')
            end 
            xline(qrs1/64); hold on; xline(T(i,1)/1024,'--g'); hold on; 
            xline(fine_T(i,1)/1024,'--b'); hold on;
            xline((AO(i,1)+window_SCG)./64); hold on 
            xline(QS2maxlined/64,'--m'); hold on;
            xline((qrs1+0.005*64)/64,'r'); hold on;
            plot(IVC(i,1)/64,IVC(i,2),'*r'); hold on
            plot(AO(i,1)/64,AO(i,2),'*r'); hold on;
            plot(RE(i,1)/64,RE(i,2),'*r'); hold on;
            plot(AC(i,1)/64,AC(i,2),'*r'); hold on;
            plot(MC(i,1)/64,MC(i,2),'*r'); hold on;
            line([AO(i,1)/64 IVC(i,1)/64],[AO(i,2) IVC(i,2)],'Color','red','LineStyle','--'); hold on;
            plot(minAO_RE(i,1)/64,minAO_RE(i,2),'*r'); hold on;
            plot(minbeforeAC(i,1)/64,minbeforeAC(i,2),'*r');
            text(IVC(i,1)/64,IVC(i,2),' IVC')
            text(AO(i,1)/64,AO(i,2),' AO')
            text(RE(i,1)/64,RE(i,2),' RE')
            text(AC(i,1)/64,AC(i,2),' AC')
            text(MC(i,1)/64,MC(i,2),' MC')
            text(minAO_RE(i,1)/64,minAO_RE(i,2),' min AO-RE');
            text(minbeforeAC(i,1)/64,minbeforeAC(i,2), 'min before AC')
            sgtitle(i)
            R_MC1(i,1)
             pause
             close all
% % %     
    %% Estraggo i parameters 
            [tIVCAO,tIVCAC,ampIVCAO,ampIVCAC,slopeIVCAO,lvet,qs2,qt,qtc,tIVCMC,tIVCRE,tIVCminAORE,tIVCminAC,ampIVCMC,ampIVCRE,ampIVCminAORE,...
                ampIVCminAC,slopeminAORERE,slopeminACAC,RdivT] = extractfeatures(64,1024,AO(i,1),AO(i,2),IVC(i,1),IVC(i,2),AC(i,1),AC(i,2),...
                RR(i-1,1),Q_SCG,fine_T(i,1),T(i,2),Q(i,1),MC(i,1),MC(i,2),RE(i,1),RE(i,2),minAO_RE(i,1),minAO_RE(i,2),minbeforeAC(i,1),...
                minbeforeAC(i,2),R(i,2));
            t_IVCAO(i,1) = tIVCAO; t_IVCAC(i,1) = tIVCAC; t_IVCRE(i,1) = tIVCRE; t_IVCminAORE(i,1) = tIVCminAORE; t_IVCminAC(i,1) = tIVCminAC;
            amp_IVCAO(i,1) = ampIVCAO; amp_IVCAC(i,1) = ampIVCAC; amp_IVCMC(i,1) = ampIVCMC; amp_IVCRE(i,1)= ampIVCRE; 
            amp_IVCminAORE(i,1) = ampIVCminAORE; amp_IVCminAC(i,1) = ampIVCminAC; R_div_T(i,1) = RdivT;
            LVET(i,1) = lvet; QS2(i,1) = qs2; QT(i,1) = qt; QTc(i,1) = qtc;
            if isnan(slopeIVCAO)
                slope_IVCAO(i,1) = 0;
            else
                slope_IVCAO(i,1) = slopeIVCAO;
            end  
            if isnan(slopeminAORERE)
                slope_minAORERE(i,1) = 0;
            else
                slope_minAORERE(i,1) = slopeminAORERE;
            end  
            if isnan(slopeminACAC)
                slope_minACAC(i,1) = 0;
            else
                slope_minACAC(i,1) = slopeminACAC;
            end  

            % QUI FACCIO IL CONTROLLO E VEDO SE DEVO METTERE ALTRI TAG O LASCIARE 0
            R_SCG = (R(i,1)/1024)*64;
            R_MC1(i,1) = R_SCG-MC(i,1); % se positivo R è dopo MC, se negativo R è prima di MC --> la maggior parte delle volte R dovrebbe essere dopo MC, quindi +
            Q_AC1(i,1) = AC(i,1)-Q_SCG;
            R_AO(i,1) = AO(i,1)-R_SCG;
    
            if R_MC1(i,1) < 0 % provo a mettere che può essere positivo
%             ma non più di 5 ms dopo R if R_MC1(i,1) < -0.005*fs_SCG
                % no sistole --> no diastole --> elimino tutto
                tag4 = tag4+1;
                R(i,3) = 4;
                MC(i,1:2) = NaN; RE(i,1:2) = NaN; AO(i,1:2) = NaN; AC(i,1:2) = NaN; IVC(i,1:2) = NaN; minAO_RE(i,1:2) = NaN; minbeforeAC(i,1:2) = NaN;
                Q(i,1:2) = NaN; fine_T(i,1:2) = NaN; T(i,1:2) = NaN; RR(i,1:2) = NaN;
                t_IVCAO(i,1) = NaN; t_IVCAC(i,1) = NaN; t_IVCMC(i,1) = NaN; t_IVCRE(i,1) = NaN; t_IVCminAORE(i,1) = NaN; t_IVCminAC(i,1) = NaN;
                amp_IVCAO(i,1) = NaN; amp_IVCAC(i,1) = NaN; amp_IVCMC(i,1) = NaN; amp_IVCRE(i,1) = NaN; amp_IVCminAORE(i,1) = NaN; amp_IVCminAC(i,1) = NaN; 
                slope_IVCAO(i,1) = NaN; slope_minAORERE(i,1) = NaN; slope_minACAC(i,1) = NaN; LVET(i,1) = NaN; QS2(i,1) = NaN; QT(i,1) = NaN; QTc(i,1) = NaN;
                LVET(i,1) = NaN; R_div_T(i,1) = NaN; R_AO(i,1) = NaN; 
            else % MC va bene ma non va bene AC
%                 if AC(i,1) > QS2maxline
                  if (Q_AC1(i,1)/64) > QS2maxline % La media sarebbe 0.44
                      % va bene la sistole ma non va bene la diastole!
                      % Cancello tutto quello che riguarda la diastole, AC
                      % e min before AC
                    tag5 = tag5+1;
                    R(i,3) = 5;
    %                 MC(i,1:2) = NaN; RE(i,1:2) = NaN; AO(i,1:2) = NaN; IVC(i,1:2) = NaN; minAO_RE(i,1:2) = NaN;
    %                 t_IVCAO(i,1) = NaN; t_IVCMC(i,1) = NaN; t_IVCRE(i,1) = NaN; t_IVCminAORE(i,1) = NaN; amp_IVCAO(i,1) = NaN; amp_IVCMC(i,1) = NaN;
%                     amp_IVCRE(i,1) = NaN; amp_IVCminAORE(i,1) = NaN; slope_IVCAO(i,1) = NaN;slope_minAORERE(i,1) = NaN;
%                     R_MC1(i,1) = NaN; QT(i,1) = NaN; QTc(i,1) = NaN;
%                     R_div_T(i,1) = NaN; R_AO(i,1) = NaN; 
                    AC(i,1:2) = NaN; minbeforeAC(i,1:2) = NaN; amp_IVCAC(i,1) = NaN;  amp_IVCminAC(i,1) = NaN; slope_minACAC(i,1) = NaN;
                    LVET(i,1) = NaN; QS2(i,1) = NaN; t_IVCAC(i,1) = NaN; t_IVCminAC(i,1) = NaN; 
                 end 
            end 

           else 
                tag2 = tag2+1;
                R(i,3) = 2;
                MC(i,1:2) = NaN; RE(i,1:2) = NaN; AO(i,1:2) = NaN; AC(i,1:2) = NaN; IVC(i,1:2) = NaN; minAO_RE(i,1:2) = NaN; minbeforeAC(i,1:2) = NaN;
                Q(i,1:2) = NaN; fine_T(i,1:2) = NaN; T(i,1:2) = NaN; RR(i,1:2) = NaN;
                t_IVCAO(i,1) = NaN; t_IVCAC(i,1) = NaN; t_IVCMC(i,1) = NaN; t_IVCRE(i,1) = NaN; t_IVCminAORE(i,1) = NaN; t_IVCminAC(i,1) = NaN;
                amp_IVCAO(i,1) = NaN; amp_IVCAC(i,1) = NaN; amp_IVCMC(i,1) = NaN; amp_IVCRE(i,1) = NaN; amp_IVCminAORE(i,1) = NaN; amp_IVCminAC(i,1) = NaN; 
                slope_IVCAO(i,1) = NaN; slope_minAORERE(i,1) = NaN; slope_minACAC(i,1) = NaN; LVET(i,1) = NaN; QS2(i,1) = NaN; QT(i,1) = NaN; QTc(i,1) = NaN;
                R_MC1(i,1) = NaN; R_AC1(i,1) = NaN; R_div_T(i,1) = NaN; R_AO(i,1) = NaN; 

        end % chiude il check dei 3 picchi dopo R
    end % chiude il numero di picchi = 1 o 2

   
    
    picchi_totali(i,1) = n_picchi(i,1);
    
    
    clearvars pks locs pksNeg locsNeg locsafterpicco locsafterR pksafterR maxpks maxlocs maxp maxsort ...
        poslocsRE locsdopoRprimaQT_ECG pksdopoRprimaQT_ECG numero locsmindopoT pksmindopoT ...
        locsmaxdopofineT pksmaxdopofineT locsprimadifineT possibileAC locsbeforeAO ...
        peakIVC locs_negdopoR pks_negdopoR locsafterIVC locsbeforeIVC pksmindopoT ...
        locsmindopoT peakminAO_RE locsbeforeAC locsprimaR_ECG pksprimaR_ECG ...
        locsdopoR_ECG pksdopoR_ECG finestraxm finestraxr p QS2maxline
end % chiude il numero dei picchi totali 

    %Tolgo le colonne che non avevo messo a zero a caso
%     FINESTREBATTITO_SCG = FINESTREBATTITO_SCG(:,1:maxlength_SCG);
%     FINESTREBATTITO_ECG = FINESTREBATTITO_ECG(:,1:maxlength_ECG);

    % Prima di salvare tolgo tutte le righe prima di iniziopicchi!!!
    MC(1:iniziopicchi-1,:) = []; RE(1:iniziopicchi-1,:) = []; AO(1:iniziopicchi-1,:) = []; AC(1:iniziopicchi-1,:) = []; R(1:iniziopicchi-1,:) = [];
    IVC(1:iniziopicchi-1,:) = []; minAO_RE(1:iniziopicchi-1,:) = []; minbeforeAC(1:iniziopicchi-1,:) = []; Q(1:iniziopicchi-1,:) = [];
    fine_T(1:iniziopicchi-1,:) = []; T(1:iniziopicchi-1,:) = []; RR(1:iniziopicchi-1,:) = [];
    t_IVCAO(1:iniziopicchi-1,:) = []; t_IVCAC(1:iniziopicchi-1,:) = []; t_IVCMC(1:iniziopicchi-1,:) = []; t_IVCRE(1:iniziopicchi-1,:) = [];
    t_IVCminAORE(1:iniziopicchi-1,:) = []; t_IVCminAC(1:iniziopicchi-1,:) = []; 
    amp_IVCAO(1:iniziopicchi-1,:) = []; amp_IVCAC(1:iniziopicchi-1,:) = []; amp_IVCMC(1:iniziopicchi-1,:) = []; amp_IVCRE(1:iniziopicchi-1,:) = []; 
    amp_IVCminAORE(1:iniziopicchi-1,:) = []; amp_IVCminAC(1:iniziopicchi-1,:) = []; R_div_T(1:iniziopicchi-1,:) = [];
    slope_IVCAO(1:iniziopicchi-1,:) = []; slope_minAORERE(1:iniziopicchi-1,:) = []; slope_minACAC(1:iniziopicchi-1,:) = [];
    LVET(1:iniziopicchi-1,:) = []; QS2(1:iniziopicchi-1,:) = []; QT(1:iniziopicchi-1,:) = []; QTc(1:iniziopicchi-1,:) = [];
%     FINESTREBATTITO_ECG(count+1:end,:) = []; FINESTREBATTITO_SCG(count+1:end,:) = []; picchi_totali(1:iniziopicchi-1,:) = [];
    R_AC1(1:iniziopicchi-1,:) = []; R_MC1(1:iniziopicchi-1,:) = []; GN(1:iniziopicchi-1,:) = []; n_picchi(1:iniziopicchi-1,:) = [];
    R_AO(1:iniziopicchi-1,:) = [];
    R_ACsec = R_AC1./64;
    R_MCsec = R_MC1./64;
    R_AOsec = R_AO./64;

    i = 0;
    for i = 1:length(R)
         tag = R(i,3);
         MC(i,3) = tag; RE(i,3) = tag; AO(i,3) = tag; AC(i,3) = tag; IVC(i,3) = tag; minAO_RE(i,3) = tag; minbeforeAC(i,3) = tag; 
         Q(i,3) = tag; fine_T(i,3) = tag; T(i,3) = tag; RR(i,3) = tag;
         t_IVCAO(i,2) = tag; t_IVCAC(i,2) = tag; t_IVCMC(i,2) = tag; t_IVCRE(i,2) = tag; t_IVCminAORE(i,2) = tag; t_IVCminAC(i,2) = tag;
         amp_IVCAO(i,2) = tag; amp_IVCAC(i,2) = tag; amp_IVCMC(i,2) = tag; amp_IVCRE(i,2) = tag; amp_IVCminAORE(i,2) = tag; 
         amp_IVCminAC(i,2) = tag; R_div_T(i,2) = tag;
         slope_IVCAO(i,2) = tag; slope_minAORERE(i,2) = tag; slope_minACAC(i,2) = tag;
         LVET(i,2) = tag; QS2(i,2) = tag; QT(i,2) = tag; QTc(i,2) = tag; R_MC(i,2) = tag; R_AC(i,2) = tag; R_div_T(i,2) = tag;
         R_AO(i,2) = tag;
    end 
        
    tag0 = length(R)-tag1-tag2-tag3-tag4-tag5-tag6-tag7-tag8;
    Perc_analizzati = tag0/length(R);
    


%     Salvo i dati (fiducial points e parameters) SU FINESTRA NORMALE
    name = erase(name,"ECG_FILT-")
    save(['C:\Users\feder\Desktop\Tesi\Data\Parameters SCG\' 'Parameters SCG-' name], 't_IVCAO','t_IVCAC','t_IVCMC','t_IVCRE','t_IVCminAORE',...
        't_IVCminAC','amp_IVCAO','amp_IVCAC','amp_IVCMC','amp_IVCRE','amp_IVCminAORE','amp_IVCminAC','R_div_T','slope_IVCAO','slope_minAORERE',...
        'slope_minACAC','LVET','QS2','QT','QTc','RR','GN','R_AO')
    save(['C:\Users\feder\Desktop\Tesi\Data\Fiducial Points SCG\' 'Fiducials SCG-' name],'MC','RE','AO','AC','R','IVC','minAO_RE','minbeforeAC',...
    'Q','fine_T','T')
%     save(['C:\Users\feder\Desktop\Tesi\Data\Windows\' 'Windows-' name],'picchi_totali','picchi_parziali','FINESTREBATTITO_ECG','FINESTREBATTITO_SCG',...
%     'R','tag0','tag1','tag2','tag3','tag4','tag5','tag6','tag7','tag8','Perc_analizzati')
    save(['C:\Users\feder\Desktop\Tesi\Data\Windows\' 'Windows-' name],'picchi_totali','picchi_parziali',...
    'R','tag0','tag1','tag2','tag3','tag4','tag5','tag6','tag7','tag8','Perc_analizzati','n_picchi')

%      Salvo i dati SU FINESTRA DI 10 SEC
%         save(['C:\Users\feder\Desktop\Tesi\Data\Parameters SCG_10SEC\' 'Parameters SCG-' name], 't_IVCAO','t_IVCAC','t_IVCMC','t_IVCRE','t_IVCminAORE',...
%         't_IVCminAC','amp_IVCAO','amp_IVCAC','amp_IVCMC','amp_IVCRE','amp_IVCminAORE','amp_IVCminAC','R_div_T','slope_IVCAO','slope_minAORERE',...
%         'slope_minACAC','LVET','QS2','QT','QTc','RR','GN','R_AO')
%     save(['C:\Users\feder\Desktop\Tesi\Data\Fiducial Points SCG_10SEC\' 'Fiducials SCG-' name],'MC','RE','AO','AC','R','IVC','minAO_RE','minbeforeAC',...
%     'Q','fine_T','T')
% %     save(['C:\Users\feder\Desktop\Tesi\Data\Windows\' 'Windows-' name],'picchi_totali','picchi_parziali','FINESTREBATTITO_ECG','FINESTREBATTITO_SCG',...
% %     'R','tag0','tag1','tag2','tag3','tag4','tag5','tag6','tag7','tag8','Perc_analizzati')
%     save(['C:\Users\feder\Desktop\Tesi\Data\Windows_10SEC\' 'Windows-' name],'picchi_totali','picchi_parziali',...
%     'R','tag0','tag1','tag2','tag3','tag4','tag5','tag6','tag7','tag8','Perc_analizzati','n_picchi')

%% Variable 'FINESTREBATTITO_ECG' was not saved. For variables larger than 2GB use MAT-file version 7.3 or later --> PROBLEMA NEL SALVATAGGIO DI FINESTRANATTITO_ECG  
end % chiude il numero di soggetti


