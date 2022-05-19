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

 end 
