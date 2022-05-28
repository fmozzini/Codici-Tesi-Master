%% Calcolo RR, AO-AO e faccio i boxplot dei vari parametri
% carico per ogni soggetto i Fiducial Points ed i Parameters 
clc 
close all
clear all

%%
folderPAR = 'C:\Users\feder\Desktop\Tesi\Data\Parameters SCG';
folderFP = 'C:\Users\feder\Desktop\Tesi\Data\Fiducial Points SCG'; 

listPAR = dir(folderPAR);
listPAR(1) = [];
listPAR(1) = [];
listFP = dir(folderFP);
listFP(1) = [];
listFP(1) = [];

 addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Parameters SCG'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Fiducial Points SCG'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Codes\
 %%
 for m = 1:1
    FOLDERPAR = fullfile(listPAR(m).folder, listPAR(m).name)
    file = dir(FOLDERPAR);
    name = file.name;
    load(name)

    FOLDERFP = fullfile(listFP(m).folder, listFP(m).name)
    file = dir(FOLDERFP);
    name = file.name;
    load(name)

    %% Calcolo RR ed AO-AO
    % RR
%     RR = diff(R(:,1)); %indici, non s o ms (divido per 1024)
%     RR_sec = RR./1024; % lo plotto in sec sulle y, 84163
% %     RR_sec(end) = [];
%     out_RR = isoutlier(RR_sec); % 84163
% %     for i = length(out_RR):-1:1
% %         if out_RR(i) == 1
% %             RR_sec(i) = [];
% %         end
% %     end 
%     lunghezza_x = length(RR_sec) 
%     x = (1:lunghezza_x)';
%     plot(x,RR_sec),xlabel('[s]'),ylabel('RR[s]'),title('RR interval')
%     boxplot(RR_sec),title('RR interval')

    %% AO-AO posizione,ampiezza
    a_count = 0;
    for a = 1:length(AO)
        if AO(a,3) == 0
            a_count = a_count + 1;
            AO_stud(a_count,:) = AO(a,:);
            R_stud(a_count,:) = R(a,:);
            amp_IVCAO_stud(a_count,:) = amp_IVCAO(a,1);
            amp_IVCAC_stud(a_count,:) = amp_IVCAC(a,1);
            t_IVCAO_stud(a_count,:) = t_IVCAO(a,1);
            t_IVCAC_stud(a_count,:) = t_IVCAC(a,1);
            slope_IVCAO_stud(a_count,:) = slope_IVCAO(a,1);
            LVET_stud(a_count,:) = LVET(a,1);
            QS2_stud(a_count,:) = QS2(a,1);
            GN_stud(a_count,:) = GN(a,1);
        end 
    end 

    q_count = 0;
    for q = 1:length(QT)
        if (AO(q,3) == 0)
            q_count = q_count + 1;
            QT_stud(q_count,:) = QT(q,1);
            QTc_stud(q_count,:) = QTc(q,1);
        end 
    end 
      
    % 
    AO_stud(end,:) = [];
    R_stud(end,:) = [];
    amp_IVCAO_stud(end,:) = [];
    amp_IVCAC_stud(end,:) = [];
    t_IVCAO_stud(end,:) = [];
    t_IVCAC_stud(end,:) = [];
    slope_IVCAO_stud(end,:) = [];
    LVET_stud(end,:) = [];
    QS2_stud(end,:) = [];
    QT_stud(end,:) = [];
    QTc_stud(end,:) = [];

    AO_AO = diff(AO_stud(:,1)); %indici, non s o ms(divido per 64) -> lungo 84163
    AO_AO_sec = AO_AO./64;
    AO_AO_tempo = AO_stud(1:end-1,1);
    AO_AO_tempo = AO_AO_tempo./64;
%     plot(AO_AO_tempo,AO_AO_sec)

%     %AO_AO - conto gli outliers
%     out = isoutlier(AO_AO_sec) % vettore con tutti gli outliers.... vale la pena andarli ad eliminare...84163
%     out = double(out); % 84163
%     % AO_AO - conto gli zeri 
%     count = 0;
%     for i = 1:length(AO_AO_sec)
%         if AO_AO_sec(i,:) == 0
%             count = count+1;
%             zeri(count) = i; %posizione degli zeri
%         end
%     end
%     zeri = zeri';
%     ZERI = zeros(length(AO_AO_sec),1); %84163
%     for i = 1:length(zeri)
%         ZERI(zeri(i),1) = 1; %VETTORE CON TUTTI 1 DOVE HO ZERO 
%     end
%     for i = 1:length(ZERI) % metto insieme una variabile in cui ho 1 se ho outliers o 0
%         if out(i,1) == 1
%             ZERI(i,1) = 1;
%         end
%     end 
%     % elimino dove ho outlier o 0 
%     for i = length(ZERI):-1:1
%         if ZERI(i) == 1
%             AO_AO_sec(i) = [];
%             AO_AO_tempo(i) = [];
%         end
%     end 

    plot(AO_AO_tempo,AO_AO_sec),xlabel('[s]'),ylabel('AO-AO[s]'),title('AO-AO interval')
    boxplot(AO_AO_sec)
%     lunghezza_x_AO_AO = length(AO_AO_sec); %70761
%     x_AO_AO = (1:lunghezza_x_AO_AO)';
%     plot(x_AO_AO,AO_AO_sec),xlabel('[s]'),ylabel('AO-AO[s]'),title('AO-AO interval')
%     boxplot(AO_AO_sec),title('AO-AO interval')

    %%
    RR = diff(R_stud(:,1)); %indici, non s o ms (divido per 1024)
    RR_sec = RR./1024; % lo plotto in sec sulle y, 84163
    RR_tempo = R_stud(1:end-1,1);
    RR_tempo = RR_tempo./1024;

%     RR_sec(end) = [];
%     out_RR = isoutlier(RR_sec); % 84163
%     for i = length(out_RR):-1:1
%         if out_RR(i) == 1
%             RR_sec(i) = [];
%         end
%     end 
%%
%     for i = length(ZERI):-1:1
%         if ZERI(i) == 1
%             RR_sec(i) = [];
%             RR_tempo(i) = [];
%         end
%     end 

    plot(RR_tempo,RR_sec),xlabel('[s]'),ylabel('RR[s]'),title('RR interval')
    boxplot(RR_sec),title('RR interval')

%     lunghezza_x = length(RR_sec) 
%     x = (1:lunghezza_x)';
%     plot(x,RR_sec),xlabel('[s]'),ylabel('RR[s]'),title('RR interval')
%     boxplot(RR_sec),title('RR interval')

    %% 
    figure()
    subplot(221), plot(RR_tempo,RR_sec),xlabel('[s]'),ylabel('R-R[s]'),title('R-R interval')
    subplot(222), plot(AO_AO_tempo,AO_AO_sec),xlabel('[s]'),ylabel('AO-AO[s]'),title('AO-AO interval')
    subplot(223), boxplot(RR_sec),title('R-R interval'),ylabel('R-R[s]')
    subplot(224), boxplot(AO_AO_sec),title('AO-AO interval'),ylabel('AO-AO[s]')
    %% PLOTTIAMO I PARAMETRI 
%     for i = length(amp_IVCAC):-1:1
%         if isnan(amp_IVCAC(i))
%             amp_IVCAC(i) = [];
%             amp_IVCAO(i) = [];
%             t_IVCAC(i) = [];
%             t_IVCAO(i) = [];
%             slope_IVCAO(i) = [];
%         end 
%     end 
       
   %%
    figure()
    subplot(221), boxplot(amp_IVCAO_stud), ylabel(['mV'])
    subplot(222), histogram(amp_IVCAO_stud), ylabel(['mV'])
    subplot(223), boxplot(amp_IVCAC_stud), ylabel(['mV'])
    subplot(224), histogram(amp_IVCAC_stud), ylabel(['mV'])

    figure()
    subplot(221), boxplot(t_IVCAO_stud), ylabel('[s]')
    subplot(222), histogram(t_IVCAO_stud), ylabel('[s]')
    subplot(223), boxplot(t_IVCAC_stud), ylabel('[s]')
    subplot(224), histogram(t_IVCAC_stud), ylabel('[s]')

    figure()
    subplot(121), boxplot(slope_IVCAO_stud)
    subplot(122), histogram(slope_IVCAO_stud)

    figure()
    subplot(221), boxplot(LVET_stud), ylabel('[s]')
    subplot(222), histogram(LVET_stud), ylabel('[s]')
    subplot(223), boxplot(QS2_stud), ylabel('[s]')
    subplot(224), histogram(QS2_stud), ylabel('[s]')

    
    notte = 0;
    giorno = 0;
    for g = 1:length(GN)
        if GN(g) == 1 %NOTTE
            notte = notte+1;
            AO_studioN(notte,:) = AO_stud(g,:);
            R_studN(notte,:) = R_stud(g,:);
            amp_IVCAO_studN(notte,:) = amp_IVCAO_stud(g,:);
            amp_IVCAC_studN(notte,:) = amp_IVCAC_stud(g,:);
            t_IVCAO_studN(notte,:) = t_IVCAO_stud(g,:);
            t_IVCAC_studN(notte,:) = t_IVCAC_stud(g,:);
            slope_IVCAO_studN(notte,:) = slope_IVCAO_stud(g,:);
            LVET_studN(notte,:) = LVET_stud(g,:);
            QS2_studN(notte,:) = QS2_stud(g,:);
            QT_studN(notte,:) = QT_stud(g,:);
            QTc_studN(notte,:) = QTc_stud(g,:);
        else % GN(i) == 0 %giorno
            giorno = giorno+1;
            AO_studioG(giorno,:) = AO_stud(g,:);
            R_studG(giorno,:) = R_stud(g,:);
            amp_IVCAO_studG(giorno,:) = amp_IVCAO_stud(g,:);
            amp_IVCAC_studG(giorno,:) = amp_IVCAC_stud(g,:);
            t_IVCAO_studG(giorno,:) = t_IVCAO_stud(g,:);
            t_IVCAC_studG(giorno,:) = t_IVCAC_stud(g,:);
            slope_IVCAO_studG(giorno,:) = slope_IVCAO_stud(g,:);
            LVET_studG(giorno,:) = LVET_stud(g,:);
            QS2_studG(giorno,:) = QS2_stud(g,:);
            QT_studG(giorno,:) = QT_stud(g,:);
            QTc_studG(giorno,:) = QTc_stud(g,:);
        end 
    end

 end 
%%
 Giorno = find(GN == 0);
 Notte = find(GN ==1);
 % Giorno
 amp_IVCAC_G = amp_IVCAC(Giorno,:);
 amp_IVCAO_G = amp_IVCAO(Giorno,:);
 amp_IVCMC_G = amp_IVCMC(Giorno,:);
 amp_IVCRE_G = amp_IVCRE(Giorno,:);
 amp_IVCminAC_G = amp_IVCminAC(Giorno,:);
 amp_IVCminAORE_G = amp_IVCminAORE(Giorno,:);
 t_IVCAC_G = t_IVCAC(Giorno,:);
 t_IVCAO_G = t_IVCAO(Giorno,:);
 t_IVCMC_G = t_IVCMC(Giorno,:);
 t_IVCRE_G = t_IVCRE(Giorno,:);
 t_IVCminAC_G = t_IVCminAC(Giorno,:);
 t_IVCminAORE_G = t_IVCminAORE(Giorno,:);
 slope_minACAC_G = slope_minACAC(Giorno,:);
 slope_minAORERE_G = slope_minAORERE(Giorno,:);
 slope_IVCAO_G = slope_IVCAO(Giorno,:);
 LVET_G = LVET(Giorno,:);
 QT_G = QT(Giorno,:);
 QTc_G = QTc(Giorno,:);
 QS2_G = QS2(Giorno,:);

 R_giorno = R(Giorno,:);
 tagG0 = 0; tagG1 = 0; tagG2 = 0; tagG3 = 0; tagG4 = 0; tagG5 = 0;
 for i  = 1:length(R_giorno)
     if R_giorno(i,3) == 1
         tagG1 = tagG1+1;
         R_giorno1(tagG1,:) = R_giorno(i,:);
     elseif R_giorno(i,3) == 2
         tagG2 = tagG2+1;
         R_giorno2(tagG2,:) = R_giorno(i,:);
     elseif R_giorno(i,3) == 3
         tagG3 = tagG3+1;
         R_giorno3(tagG3,:) = R_giorno(i,:);
     elseif R_giorno(i,3) == 4
         tagG4 = tagG4+1;
         R_giorno4(tagG4,:) = R_giorno(i,:);
     elseif R_giorno(i,3) == 5
         tagG5 = tagG5+1;
         R_giorno5(tagG5,:) = R_giorno(i,:);
     else 
         tagG0 = tagG0+1;
         R_giorno0(tagG0,:) = R_giorno(i,:);
         amp_IVCAC_G0(tagG0,1) = amp_IVCAC_G(i,1);
         amp_IVCAO_G0(tagG0,1) = amp_IVCAO_G(i,1);
         amp_IVCMC_G0(tagG0,1) = amp_IVCMC_G(i,1);
         amp_IVCRE_G0(tagG0,1) = amp_IVCRE_G(i,1);
         amp_IVCminAC_G0(tagG0,1) = amp_IVCminAC_G(i,1);
         amp_IVCminAORE_G0(tagG0,1) = amp_IVCminAORE_G(i,1);
         t_IVCAC_G0(tagG0,1) = t_IVCAC_G(i,1);
         t_IVCAO_G0(tagG0,1) = t_IVCAO_G(i,1);
         t_IVCMC_G0(tagG0,1) = t_IVCMC_G(i,1);
         t_IVCRE_G0(tagG0,1) = t_IVCRE_G(i,1);
         t_IVCminAC_G0(tagG0,1) = t_IVCminAC_G(i,1);
         t_IVCminAORE_G0(tagG0,1) = t_IVCminAORE_G(i,1);
         slope_IVCAO_G0(tagG0,1) = slope_IVCAO_G(i,1);
         slope_minACAC_G0(tagG0,1) = slope_minACAC_G(i,1);
         slope_minAORERE_G0(tagG0,1) = slope_minAORERE_G(i,1);
         LVET_G0(tagG0,1) = LVET_G(i,1);
         QT_G0(tagG0,1) = QT_G(i,1);
         QTc_G0(tagG0,1) = QTc_G(i,1);
         QS2_G0(tagG0,1) = QS2_G(i,1);
     end 
 end 
% Notte
 amp_IVCAC_N = amp_IVCAC(Notte,:);
 amp_IVCAO_N = amp_IVCAO(Notte,:);
 amp_IVCMC_N = amp_IVCMC(Notte,:);
 amp_IVCRE_N = amp_IVCRE(Notte,:);
 amp_IVCminAC_N = amp_IVCminAC(Notte,:);
 amp_IVCminAORE_N = amp_IVCminAORE(Notte,:);
 t_IVCAC_N = t_IVCAC(Notte,:);
 t_IVCAO_N = t_IVCAO(Notte,:);
 t_IVCMC_N = t_IVCMC(Notte,:);
 t_IVCRE_N = t_IVCRE(Notte,:);
 t_IVCminAC_N = t_IVCminAC(Notte,:);
 t_IVCminAORE_N = t_IVCminAORE(Notte,:);
 slope_IVCAO_N = slope_IVCAO(Notte,:);
 slope_minACAC_N = slope_minACAC(Notte,:);
 slope_minAORERE_N = slope_minAORERE(Notte,:);
 LVET_N = LVET(Notte,:);
 QT_N = QT(Notte,:);
 QTc_N = QTc(Notte,:);
 QS2_N = QS2(Notte,:);

 R_notte = R(Notte,:);
 tagN0 = 0; tagN1 = 0; tagN2 = 0; tagN3 = 0; tagN4 = 0; tagN5 = 0;
 for i  = 1:length(R_notte)
     if R_notte(i,3) == 1
         tagN1 = tagN1+1;
         R_notte1(tagN1,:) = R_notte(i,3);
     elseif R_notte(i,3) == 2
         tagN2 = tagN2+1;
         R_notte2(tagN2,:) = R_notte(i,3);
     elseif R_notte(i,3) == 3
         tagN3 = tagN3+1;
         R_notte3(tagN3,:) = R_notte(i,3);
     elseif R_notte(i,3) == 4
         tagN4 = tagN4+1;
         R_notte4(tagN4,:) = R_notte(i,3);
     elseif R_notte(i,3) == 5
         tagN5 = tagN5+1;
         R_notte5(tagN5,:) = R_notte(i,3);
     else 
         tagN0 = tagN0+1;
         R_notte0(tagN0,:) = R_notte(i,3);
         amp_IVCAC_N0(tagN0,1) = amp_IVCAC_N(i,1);
         amp_IVCAO_N0(tagN0,1) = amp_IVCAO_N(i,1);
         amp_IVCMC_N0(tagN0,1) = amp_IVCMC_N(i,1);
         amp_IVCRE_N0(tagN0,1) = amp_IVCRE_N(i,1);
         amp_IVCminAC_N0(tagN0,1) = amp_IVCminAC_N(i,1);
         amp_IVCminAORE_N0(tagN0,1) = amp_IVCminAORE_N(i,1);
         t_IVCAC_N0(tagN0,1) = t_IVCAC_N(i,1);
         t_IVCAO_N0(tagN0,1) = t_IVCAO_N(i,1);
         t_IVCMC_N0(tagN0,1) = t_IVCMC_N(i,1);
         t_IVCRE_N0(tagN0,1) = t_IVCRE_N(i,1);
         t_IVCminAC_N0(tagN0,1) = t_IVCminAC_N(i,1);
         t_IVCminAORE_N0(tagN0,1) = t_IVCminAORE_N(i,1);
         slope_IVCAO_N0(tagN0,1) = slope_IVCAO_N(i,1);
         slope_minACAC_N0(tagN0,1) = slope_minACAC_N(i,1);
         slope_minAORERE_N0(tagN0,1) = slope_minAORERE_N(i,1);
         LVET_N0(tagN0,1) = LVET_N(i,1);
         QT_N0(tagN0,1) = QT_N(i,1);
         QTc_N0(tagN0,1) = QTc_N(i,1);
         QS2_N0(tagN0,1) = QS2_N(i,1);
     end 
 end 
%%
 % Histogram - Amplitudes 
 figure()
 subplot(321),histogram(amp_IVCAO_G0); hold on; histogram(amp_IVCAO_N0),xlabel(['mV']),title('Amplitude IVC-AO'),legend('Day','Night')
 subplot(322),histogram(amp_IVCAC_G0); hold on; histogram(amp_IVCAC_N0),xlabel(['mV']),title('Amplitude IVC-AC'),legend('Day','Night')
 subplot(323),histogram(amp_IVCMC_G0); hold on; histogram(amp_IVCMC_N0),xlabel(['mV']),title('Amplitude IVC-MC'),legend('Day','Night')
 subplot(324),histogram(amp_IVCRE_G0); hold on; histogram(amp_IVCRE_N0),xlabel(['mV']),title('Amplitude IVC-RE'),legend('Day','Night')
 subplot(325),histogram(amp_IVCminAC_G0); hold on; histogram(amp_IVCminAC_N0),xlabel(['mV']),title('Amplitude IVC-minimum before AC'),legend('Day','Night')
 subplot(326),histogram(amp_IVCminAORE_G0); hold on; histogram(amp_IVCminAORE_N0),xlabel(['mV']),title('Amplitude IVC-minimum between AO and RE'),legend('Day','Night')

%% Histogram - Time intervals 
 figure()
 subplot(321),histogram(t_IVCAO_G0./1000); hold on; histogram(t_IVCAO_N0./1000),xlabel(['ms']),title('Time interval IVC-AO'),legend('Day','Night')
 subplot(322),histogram(t_IVCAC_G0./1000); hold on; histogram(t_IVCAC_N0./1000),xlabel(['ms']),title('Time interval IVC-AC'),legend('Day','Night')
 subplot(323),histogram(t_IVCMC_G0./1000); hold on; histogram(t_IVCMC_N0./1000),xlabel(['ms']),title('Time interval IVC-MC'),legend('Day','Night')
 subplot(324),histogram(t_IVCRE_G0./1000); hold on; histogram(t_IVCRE_N0./1000),xlabel(['ms']),title('Time interval IVC-RE'),legend('Day','Night')
 subplot(325),histogram(t_IVCminAC_G0./1000); hold on; histogram(t_IVCminAC_N0./1000),xlabel(['ms']),title('Time interval IVC-minimum before AC'),legend('Day','Night')
 subplot(326),histogram(t_IVCminAORE_G0./1000); hold on; histogram(t_IVCminAORE_N0./1000),xlabel(['ms']),title('Time interval IVC-minimum between AO and RE'),legend('Day','Night')
 
 %% Histogram - Slope
 figure()
 subplot(221),histogram(slope_IVCAO_G0); hold on; histogram(slope_IVCAO_N0),title('Slope IVC-AO'),legend('Day','Night')
 subplot(222),histogram(slope_minAORERE_G0); hold on; histogram(slope_minAORERE_N0),title('Slope minimum between AO and RE-RE'),legend('Day','Night')
 subplot(223),histogram(slope_minACAC_G0); hold on; histogram(slope_minACAC_N0),title('Slope minimum before AC-AC'),legend('Day','Night')
 
 %% Histogram - Other parameters
 figure()
 subplot(221),histogram(LVET_G0./1000); hold on; histogram(LVET_N0./1000),xlabel(['ms']),title('Left Ventricular Ejection Time'),legend('Day','Night')
 subplot(222),histogram(QS2_G0./1000); hold on; histogram(QS2_N0./1000),xlabel(['ms']),title('QS2'),legend('Day','Night')
 subplot(223),histogram(QT_G0./1000); hold on; histogram(QT_N0./1000),xlabel(['ms']),title('QT'),legend('Day','Night')
 subplot(224),histogram(QTc_G0./1000); hold on; histogram(QTc_N0./1000),xlabel(['ms']),title('Corrected QT'),legend('Day','Night')
%%
% Forse dovrei aggiungere anche le info di Cinque ! BOHHHH
 figure()
    subplot(221), boxplot(amp_IVCAO_G0), ylabel(['mV'])
    subplot(222), histogram(amp_IVCAO_stud), ylabel(['mV'])
    subplot(223), boxplot(amp_IVCAC_stud), ylabel(['mV'])
    subplot(224), histogram(amp_IVCAC_stud), ylabel(['mV'])

    figure()
    subplot(221), boxplot(t_IVCAO_stud), ylabel('[s]')
    subplot(222), histogram(t_IVCAO_stud), ylabel('[s]')
    subplot(223), boxplot(t_IVCAC_stud), ylabel('[s]')
    subplot(224), histogram(t_IVCAC_stud), ylabel('[s]')

    figure()
    subplot(121), boxplot(slope_IVCAO_stud)
    subplot(122), histogram(slope_IVCAO_stud)

    figure()
    subplot(221), boxplot(LVET_stud), ylabel('[s]')
    subplot(222), histogram(LVET_stud), ylabel('[s]')
    subplot(223), boxplot(QS2_stud), ylabel('[s]')
    subplot(224), histogram(QS2_stud), ylabel('[s]')


