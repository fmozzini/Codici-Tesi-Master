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
    RR = diff(R(:,1)); %indici, non s o ms (divido per 1024)
    RR_sec = RR./1024; % lo plotto in sec sulle y, 84163
%     RR_sec(end) = [];
    out_RR = isoutlier(RR_sec); % 84163
    for i = length(out_RR):-1:1
        if out_RR(i) == 1
            RR_sec(i) = [];
        end
    end 
%     RR_sec(84074) = [];  % togliere riga 84074 - è sicuramente un outlier
    lunghezza_x = length(RR_sec) 
    x = (1:lunghezza_x)';
    plot(x,RR_sec),xlabel('[s]'),ylabel('RR[s]'),title('RR interval')


    %% AO-AO posizione,ampiezza
    AO_AO = diff(AO(:,1)); %indici, non s o ms(divido per 64) -> lungo 84163
    AO_AO_sec = AO_AO./64;

    % AO_AO - conto gli outliers
    out = isoutlier(AO_AO_sec) % vettore con tutti gli outliers.... vale la pena andarli ad eliminare...84163
    out = double(out); % 84163
    % AO_AO - conto gli zeri 
    count = 0;
    for i = 1:length(AO_AO_sec)
        if AO_AO_sec(i,:) == 0
            count = count+1;
            zeri(count) = i; %posizione degli zeri
        end
    end
    zeri = zeri';
    ZERI = zeros(length(AO_AO_sec),1); %84163
    for i = 1:length(zeri)
        ZERI(zeri(i),1) = 1; %VETTORE CON TUTTI 1 DOVE HO ZERO 
    end
    for i = 1:length(ZERI) % metto insieme una variabile in cui ho 1 se ho outliers o 0
        if out(i,1) == 1
            ZERI(i,1) = 1;
        end
    end 
    % elimino dove ho outlier o 0 
    for i = length(ZERI):-1:1
        if ZERI(i) == 1
            AO_AO_sec(i) = [];
        end
    end 

    lunghezza_x_AO_AO = length(AO_AO_sec); %70761
    x_AO_AO = (1:lunghezza_x_AO_AO)';
    plot(x_AO_AO,AO_AO_sec),xlabel('[s]'),ylabel('AO-AO[s]'),title('AO-AO interval')

    %% 
    figure()
    subplot(211), plot(x,RR_sec),xlabel('[s]'),ylabel('RR[s]'),title('RR interval')
    subplot(212), plot(x_AO_AO,AO_AO_sec),xlabel('[s]'),ylabel('AO-AO[s]'),title('AO-AO interval')
    
    %% PLOTTIAMO I PARAMETRI 
    figure()
    boxplot(amp_IVCAO),title('amp IVC-AO [??]')
    figure()
    boxplot(amp_IVCAC),title('amp IVC-AC [??]')
    figure()
    boxplot(t_IVCAO),title('t IVC-AO [s]')
    figure()
    boxplot(t_IVCAC),title('t IVC-AC [s]')
    figure()
    boxplot(slope_IVCAO),title('slope IVC-AO')

 end 
