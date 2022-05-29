%% Federica Mozzini - 946400
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
     tagG0 = 0; tagG1 = 0; tagG2 = 0; tagG3 = 0; tagG4 = 0; tagG5 = 0; tagG6 = 0;
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
         elseif R_giorno(i,3) == 6
             tagG6 = tagG6+1;
             R_giorno6(tagG6,:) = R_giorno(i,:);
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
     tagN0 = 0; tagN1 = 0; tagN2 = 0; tagN3 = 0; tagN4 = 0; tagN5 = 0; tagN6 = 0;
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
         elseif R_notte(i,3) == 6
             tagN6 = tagN6+1;
             R_notte6(tagN6,:) = R_notte(i,3);
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
     
%      Num_G05 = length(R_giorno0)+length(R_giorno5);
     tagG_0_5 = 0;
     for i = 1:length(R_giorno)
         if R_giorno(i,3) == 0 || R_giorno(i,3) == 5
             tagG_0_5 = tagG_0_5+1;
             R_giorno0_5(tagG_0_5,:) = R_giorno(i,:);
             amp_IVCAC_G_0_5(tagG_0_5,1) = amp_IVCAC_G(i,1);
             amp_IVCAO_G_0_5(tagG_0_5,1) = amp_IVCAO_G(i,1);
             amp_IVCMC_G_0_5(tagG_0_5,1) = amp_IVCMC_G(i,1);
             amp_IVCRE_G_0_5(tagG_0_5,1) = amp_IVCRE_G(i,1);
             amp_IVCminAC_G_0_5(tagG_0_5,1) = amp_IVCminAC_G(i,1);
             amp_IVCminAORE_G_0_5(tagG_0_5,1) = amp_IVCminAORE_G(i,1);
             t_IVCAC_G_0_5(tagG_0_5,1) = t_IVCAC_G(i,1);
             t_IVCAO_G_0_5(tagG_0_5,1) = t_IVCAO_G(i,1);
             t_IVCMC_G_0_5(tagG_0_5,1) = t_IVCMC_G(i,1);
             t_IVCRE_G_0_5(tagG_0_5,1) = t_IVCRE_G(i,1);
             t_IVCminAC_G_0_5(tagG_0_5,1) = t_IVCminAC_G(i,1);
             t_IVCminAORE_G_0_5(tagG_0_5,1) = t_IVCminAORE_G(i,1);
             slope_IVCAO_G_0_5(tagG_0_5,1) = slope_IVCAO_G(i,1);
             slope_minACAC_G_0_5(tagG_0_5,1) = slope_minACAC_G(i,1);
             slope_minAORERE_G_0_5(tagG_0_5,1) = slope_minAORERE_G(i,1);
             LVET_G_0_5(tagG_0_5,1) = LVET_G(i,1);
             QT_G_0_5(tagG_0_5,1) = QT_G(i,1);
             QTc_G_0_5(tagG_0_5,1) = QTc_G(i,1);
             QS2_G_0_5(tagG_0_5,1) = QS2_G(i,1);
         end 
     end 

     tagN_0_5 = 0;
     for i = 1:length(R_notte)
         if R_notte(i,3) == 0 || R_notte(i,3) == 5
             tagN_0_5 = tagN_0_5+1;
             R_notte0_5(tagN_0_5,:) = R_notte(i,:);
             amp_IVCAC_N_0_5(tagN_0_5,1) = amp_IVCAC_N(i,1);
             amp_IVCAO_N_0_5(tagN_0_5,1) = amp_IVCAO_N(i,1);
             amp_IVCMC_N_0_5(tagN_0_5,1) = amp_IVCMC_N(i,1);
             amp_IVCRE_N_0_5(tagN_0_5,1) = amp_IVCRE_N(i,1);
             amp_IVCminAC_N_0_5(tagN_0_5,1) = amp_IVCminAC_N(i,1);
             amp_IVCminAORE_N_0_5(tagN_0_5,1) = amp_IVCminAORE_N(i,1);
             t_IVCAC_N_0_5(tagN_0_5,1) = t_IVCAC_N(i,1);
             t_IVCAO_N_0_5(tagN_0_5,1) = t_IVCAO_N(i,1);
             t_IVCMC_N_0_5(tagN_0_5,1) = t_IVCMC_N(i,1);
             t_IVCRE_N_0_5(tagN_0_5,1) = t_IVCRE_N(i,1);
             t_IVCminAC_N_0_5(tagN_0_5,1) = t_IVCminAC_N(i,1);
             t_IVCminAORE_N_0_5(tagN_0_5,1) = t_IVCminAORE_N(i,1);
             slope_IVCAO_N_0_5(tagN_0_5,1) = slope_IVCAO_N(i,1);
             slope_minACAC_N_0_5(tagN_0_5,1) = slope_minACAC_N(i,1);
             slope_minAORERE_N_0_5(tagN_0_5,1) = slope_minAORERE_N(i,1);
             LVET_N_0_5(tagN_0_5,1) = LVET_N(i,1);
             QT_N_0_5(tagN_0_5,1) = QT_N(i,1);
             QTc_N_0_5(tagN_0_5,1) = QTc_N(i,1);
             QS2_N_0_5(tagN_0_5,1) = QS2_N(i,1);
         end 
     end 
     % Faccio un nuovo caso in cui ho sia il tag0 che il tag5
%%
     % Histogram - Amplitudes - TAG 0
     figure()
     subplot(321),histogram(amp_IVCAO_G0); hold on; histogram(amp_IVCAO_N0),xlabel(['mV']),title('Amplitude IVC-AO'),legend('Day','Night')
     subplot(322),histogram(amp_IVCAC_G0); hold on; histogram(amp_IVCAC_N0),xlabel(['mV']),title('Amplitude IVC-AC'),legend('Day','Night')
     subplot(323),histogram(amp_IVCMC_G0); hold on; histogram(amp_IVCMC_N0),xlabel(['mV']),title('Amplitude IVC-MC'),legend('Day','Night')
     subplot(324),histogram(amp_IVCRE_G0); hold on; histogram(amp_IVCRE_N0),xlabel(['mV']),title('Amplitude IVC-RE'),legend('Day','Night')
     subplot(325),histogram(amp_IVCminAC_G0); hold on; histogram(amp_IVCminAC_N0),xlabel(['mV']),title('Amplitude IVC-minimum before AC'),legend('Day','Night')
     subplot(326),histogram(amp_IVCminAORE_G0); hold on; histogram(amp_IVCminAORE_N0),xlabel(['mV']),title('Amplitude IVC-minimum between AO and RE'),legend('Day','Night')
     sgtitle('Amplitudes - Tag 0')

    % Histogram - Time intervals - TAG 0
     figure()
     subplot(321),histogram(t_IVCAO_G0./1000); hold on; histogram(t_IVCAO_N0./1000),xlabel(['ms']),title('Time interval IVC-AO'),legend('Day','Night')
     subplot(322),histogram(t_IVCAC_G0./1000); hold on; histogram(t_IVCAC_N0./1000),xlabel(['ms']),title('Time interval IVC-AC'),legend('Day','Night')
     subplot(323),histogram(t_IVCMC_G0./1000); hold on; histogram(t_IVCMC_N0./1000),xlabel(['ms']),title('Time interval IVC-MC'),legend('Day','Night')
     subplot(324),histogram(t_IVCRE_G0./1000); hold on; histogram(t_IVCRE_N0./1000),xlabel(['ms']),title('Time interval IVC-RE'),legend('Day','Night')
     subplot(325),histogram(t_IVCminAC_G0./1000); hold on; histogram(t_IVCminAC_N0./1000),xlabel(['ms']),title('Time interval IVC-minimum before AC'),legend('Day','Night')
     subplot(326),histogram(t_IVCminAORE_G0./1000); hold on; histogram(t_IVCminAORE_N0./1000),xlabel(['ms']),title('Time interval IVC-minimum between AO and RE'),legend('Day','Night')
     sgtitle('Time Intervals - Tag 0')
 
     % Histogram - Slopes - TAG 0
     figure()
     subplot(221),histogram(slope_IVCAO_G0); hold on; histogram(slope_IVCAO_N0),title('Slope IVC-AO'),legend('Day','Night')
     subplot(222),histogram(slope_minAORERE_G0); hold on; histogram(slope_minAORERE_N0),title('Slope minimum between AO and RE-RE'),legend('Day','Night')
     subplot(223),histogram(slope_minACAC_G0); hold on; histogram(slope_minACAC_N0),title('Slope minimum before AC-AC'),legend('Day','Night')
     sgtitle('Slopes - Tag 0')
 
     % Histogram - Other parameters - TAG 0
     figure()
     subplot(221),histogram(LVET_G0./1000); hold on; histogram(LVET_N0./1000),xlabel(['ms']),title('Left Ventricular Ejection Time'),legend('Day','Night')
     subplot(222),histogram(QS2_G0./1000); hold on; histogram(QS2_N0./1000),xlabel(['ms']),title('QS2'),legend('Day','Night')
     subplot(223),histogram(QT_G0./1000); hold on; histogram(QT_N0./1000),xlabel(['ms']),title('QT'),legend('Day','Night')
     subplot(224),histogram(QTc_G0./1000); hold on; histogram(QTc_N0./1000),xlabel(['ms']),title('Corrected QT'),legend('Day','Night')
     sgtitle('Other Parameters - Tag 0')

      %% Histogram - Amplitudes - TAG 0 E 5
     figure()
     subplot(321),histogram(amp_IVCAO_G_0_5); hold on; histogram(amp_IVCAO_N_0_5),xlabel(['mV']),title('Amplitude IVC-AO'),legend('Day','Night')
     subplot(322),histogram(amp_IVCAC_G_0_5); hold on; histogram(amp_IVCAC_N_0_5),xlabel(['mV']),title('Amplitude IVC-AC'),legend('Day','Night')
     subplot(323),histogram(amp_IVCMC_G_0_5); hold on; histogram(amp_IVCMC_N_0_5),xlabel(['mV']),title('Amplitude IVC-MC'),legend('Day','Night')
     subplot(324),histogram(amp_IVCRE_G_0_5); hold on; histogram(amp_IVCRE_N_0_5),xlabel(['mV']),title('Amplitude IVC-RE'),legend('Day','Night')
     subplot(325),histogram(amp_IVCminAC_G_0_5); hold on; histogram(amp_IVCminAC_N_0_5),xlabel(['mV']),title('Amplitude IVC-minimum before AC'),legend('Day','Night')
     subplot(326),histogram(amp_IVCminAORE_G_0_5); hold on; histogram(amp_IVCminAORE_N_0_5),xlabel(['mV']),title('Amplitude IVC-minimum between AO and RE'),legend('Day','Night')
     sgtitle('Amplitudes - Tag 0 & 5')

     % Histogram - Time intervals  - TAG 0 E 5 
     figure()
     subplot(321),histogram(t_IVCAO_G_0_5./1000); hold on; histogram(t_IVCAO_N_0_5./1000),xlabel(['ms']),title('Time interval IVC-AO'),legend('Day','Night')
     subplot(322),histogram(t_IVCAC_G_0_5./1000); hold on; histogram(t_IVCAC_N_0_5./1000),xlabel(['ms']),title('Time interval IVC-AC'),legend('Day','Night')
     subplot(323),histogram(t_IVCMC_G_0_5./1000); hold on; histogram(t_IVCMC_N_0_5./1000),xlabel(['ms']),title('Time interval IVC-MC'),legend('Day','Night')
     subplot(324),histogram(t_IVCRE_G_0_5./1000); hold on; histogram(t_IVCRE_N_0_5./1000),xlabel(['ms']),title('Time interval IVC-RE'),legend('Day','Night')
     subplot(325),histogram(t_IVCminAC_G_0_5./1000); hold on; histogram(t_IVCminAC_N_0_5./1000),xlabel(['ms']),title('Time interval IVC-minimum before AC'),legend('Day','Night')
     subplot(326),histogram(t_IVCminAORE_G_0_5./1000); hold on; histogram(t_IVCminAORE_N_0_5./1000),xlabel(['ms']),title('Time interval IVC-minimum between AO and RE'),legend('Day','Night')
     sgtitle('Time Intervals - Tag 0 & 5')
 
     % Histogram - Slope - TAG 0 E 5
     figure()
     subplot(221),histogram(slope_IVCAO_G_0_5); hold on; histogram(slope_IVCAO_N_0_5),title('Slope IVC-AO'),legend('Day','Night')
     subplot(222),histogram(slope_minAORERE_G_0_5); hold on; histogram(slope_minAORERE_N_0_5),title('Slope minimum between AO and RE-RE'),legend('Day','Night')
     subplot(223),histogram(slope_minACAC_G_0_5); hold on; histogram(slope_minACAC_N_0_5),title('Slope minimum before AC-AC'),legend('Day','Night')
     sgtitle('Slopes - Tag 0 & 5')
 
     % Histogram - Other parameters - TAG 0 E 5
     figure()
     subplot(221),histogram(LVET_G_0_5./1000); hold on; histogram(LVET_N_0_5./1000),xlabel(['ms']),title('Left Ventricular Ejection Time'),legend('Day','Night')
     subplot(222),histogram(QS2_G_0_5./1000); hold on; histogram(QS2_N_0_5./1000),xlabel(['ms']),title('QS2'),legend('Day','Night')
     subplot(223),histogram(QT_G_0_5./1000); hold on; histogram(QT_N_0_5./1000),xlabel(['ms']),title('QT'),legend('Day','Night')
     subplot(224),histogram(QTc_G_0_5./1000); hold on; histogram(QTc_N_0_5./1000),xlabel(['ms']),title('Corrected QT'),legend('Day','Night')
     sgtitle('Other Parameters - Tag 0 & 5')
 end 
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

