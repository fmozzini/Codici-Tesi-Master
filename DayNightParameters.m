%% Federica Mozzini - 946400
%% Calcolo RR, AO-AO e faccio i boxplot dei vari parametri
% carico per ogni soggetto i Fiducial Points ed i Parameters 
clc 
close all
clear all

%%
% folderPAR = 'C:\Users\feder\Desktop\Tesi\Data\Parameters SCG';
% folderFP = 'C:\Users\feder\Desktop\Tesi\Data\Fiducial Points SCG'; 
folderPAR = 'C:\Users\feder\Desktop\Tesi\Data\Parameters SCG_10SEC';
folderFP = 'C:\Users\feder\Desktop\Tesi\Data\Fiducial Points SCG_10SEC'; 

listPAR = dir(folderPAR);
listPAR(1) = [];
listPAR(1) = [];
listFP = dir(folderFP);
listFP(1) = [];
listFP(1) = [];

%  addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Parameters SCG'\
%  addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Fiducial Points SCG'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Parameters SCG_10SEC'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Fiducial Points SCG_10SEC'\
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
     PEP_G = PEP(Giorno,:);
     R_MC1_G = R_MC1(Giorno,:);
     R_AO_G = R_AO(Giorno,:);
     R_AC1_G = R_AC1(Giorno,:);
    
     R_giorno = R(Giorno,:);
     tagG0 = 0; tagG1 = 0; tagG2 = 0; tagG3 = 0; tagG4 = 0; tagG5 = 0; tagG6 = 0; tagG7 = 0; tagG8 = 0;
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
         elseif R_giorno(i,3) == 7
             tagG7 = tagG7+1;
             R_giorno7(tagG7,:) = R_giorno(i,:);
         elseif R_giorno(i,3) == 8
             tagG8 = tagG8+1;
             R_giorno8(tagG8,:) = R_giorno(i,:);
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
             PEP_G0(tagG0,1) = PEP_G(i,1);
             R_MC1_G0(tagG0,1) = R_MC1_G(i,1);
             R_AO_G0(tagG0,1) = R_AO_G(i,1);
             R_AC1_G0(tagG0,1) = R_AC1_G(i,1);
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
     PEP_N = PEP(Notte,:);
     R_AO_N = R_AO(Notte,:);
     R_AC1_N = R_AC1(Notte,:);
     R_MC1_N = R_MC1(Notte,:);
    
     R_notte = R(Notte,:);
     tagN0 = 0; tagN1 = 0; tagN2 = 0; tagN3 = 0; tagN4 = 0; tagN5 = 0; tagN6 = 0; tagN7 = 0; tagN8 = 0;
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
         elseif R_notte(i,3) == 7
             tagN7 = tagN7+1;
             R_notte7(tagN7,:) = R_notte(i,3);
         elseif R_notte(i,3) == 8
             tagN8 = tagN8+1;
             R_notte8(tagN8,:) = R_notte(i,3);
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
             PEP_N0(tagN0,1) = PEP_N(i,1);
             R_AO_N0(tagN0,1) = R_AO_N(i,1);
             R_AC1_N0(tagN0,1) = R_AC1_N(i,1);
             R_MC1_N0(tagN0,1) = R_MC1_N(i,1);
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
             PEP_G_0_5(tagG_0_5,1) = PEP_G(i,1);
             R_AO_G_0_5(tagG_0_5,1) = R_AO_G(i,1);
             R_AC1_G_0_5(tagG_0_5,1) = R_AC1_G(i,1);
             R_MC1_G_0_5(tagG_0_5,1) = R_MC1_G(i,1);
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
             PEP_N_0_5(tagN_0_5,1) = PEP_N(i,1);
             R_AO_N_0_5(tagN_0_5,1) = R_AO_N(i,1);
             R_AC1_N_0_5(tagN_0_5,1) = R_AC1_N(i,1);
             R_MC1_N_0_5(tagN_0_5,1) = R_MC1_N(i,1);
         end 
     end 

     % Aggiusto nel caso di tag0, tolto solo i minori di zero 
     minorizero = find(LVET_G0<0);
     for l = length(minorizero):-1:1
         LVET_G0(minorizero(l)) = [];
     end 
     clear minorizero
     minorizero = find(LVET_N0<0);
     for l = length(minorizero):-1:1
         LVET_N0(minorizero(l)) = [];
     end
     clear minorizero
     minorizero = find(QS2_G0<0);
     for l = length(minorizero):-1:1
         QS2_G0(minorizero(l)) = [];
     end 
     clear minorizero
     minorizero = find(QS2_N0<0);
     for l = length(minorizero):-1:1
         QS2_N0(minorizero(l)) = [];
     end 
     clear minorizero

     minorizero = find(t_IVCAC_G0<0);
     for l = length(minorizero):-1:1
         t_IVCAC_G0(minorizero(l)) = [];
     end 
     clear minorizero
     minorizero = find(t_IVCAC_N0<0);
     for l = length(minorizero):-1:1
         t_IVCAC_N0(minorizero(l)) = [];
     end
     clear minorizero

     minorizero = find(t_IVCminAC_G0<0);
     for l = length(minorizero):-1:1
         t_IVCminAC_G0(minorizero(l)) = [];
     end 
     clear minorizero
     minorizero = find(t_IVCminAC_N0<0);
     for l = length(minorizero):-1:1
         t_IVCminAC_N0(minorizero(l)) = [];
     end
     clear minorizero

      minorizero = find(PEP_G0<0);
     for l = length(minorizero):-1:1
         PEP_G0(minorizero(l)) = [];
     end 
     clear minorizero
     minorizero = find(PEP_N0<0);
     for l = length(minorizero):-1:1
         PEP_N0(minorizero(l)) = [];
     end
     clear minorizero

      minorizero = find(R_AO_G0<0);
     for l = length(minorizero):-1:1
         R_AO_G0(minorizero(l)) = [];
     end 
     clear minorizero
     minorizero = find(R_AO_N0<0);
     for l = length(minorizero):-1:1
         R_AO_N0(minorizero(l)) = [];
     end
     clear minorizero

      minorizero = find(R_AC1_G0<0);
     for l = length(minorizero):-1:1
         R_AC1_G0(minorizero(l)) = [];
     end 
     clear minorizero
     minorizero = find(R_AC1_N0<0);
     for l = length(minorizero):-1:1
         R_AC1_N0(minorizero(l)) = [];
     end
     clear minorizero

      minorizero = find(R_MC1_G0<0);
     for l = length(minorizero):-1:1
         R_MC1_G0(minorizero(l)) = [];
     end 
     clear minorizero
     minorizero = find(R_MC1_N0<0);
     for l = length(minorizero):-1:1
         R_MC1_N0(minorizero(l)) = [];
     end
     clear minorizero

     % AGGIUSTARE DA TUTTI I LVET E QS2 TOGLIENDO LE RIGHE DEL NAN
     for l = length(LVET_G_0_5):-1:1
         if isnan(LVET_G_0_5(l))
             LVET_G_0_5(l) = [];
         end
     end
     minorizero = find(LVET_G_0_5<0);
     for l = length(minorizero):-1:1
         LVET_G_0_5(minorizero(l)) = [];
     end 
     clear minorizero

     for l = length(LVET_N_0_5):-1:1
         if isnan(LVET_N_0_5(l))
             LVET_N_0_5(l) = [];
         end
     end
     minorizero = find(LVET_N_0_5<0);
     for l = length(minorizero):-1:1
         LVET_N_0_5(minorizero(l)) = [];
     end 
     clear minorizero
     %QS2
    for l = length(QS2_G_0_5):-1:1
         if isnan(QS2_G_0_5(l))
             QS2_G_0_5(l) = [];
         end
     end
     minorizero = find(QS2_G_0_5<0);
     for l = length(minorizero):-1:1
         QS2_G_0_5(minorizero(l)) = [];
     end 
     clear minorizero

     for l = length(QS2_N_0_5):-1:1
         if isnan(QS2_N_0_5(l))
             QS2_N_0_5(l) = [];
         end
     end
     minorizero = find(QS2_N_0_5<0);
     for l = length(minorizero):-1:1
         QS2_N_0_5(minorizero(l)) = [];
     end 
     clear minorizero

     % Tempo IVC-AC
     for l = length(t_IVCAC_G_0_5):-1:1
         if isnan(t_IVCAC_G_0_5(l))
             t_IVCAC_G_0_5(l) = [];
         end
     end
     minorizero = find(t_IVCAC_G_0_5<0);
     for l = length(minorizero):-1:1
         t_IVCAC_G_0_5(minorizero(l)) = [];
     end 
     clear minorizero

     for l = length(t_IVCAC_N_0_5):-1:1
         if isnan(t_IVCAC_N_0_5(l))
             t_IVCAC_N_0_5(l) = [];
         end
     end
     minorizero = find(t_IVCAC_N_0_5<0);
     for l = length(minorizero):-1:1
         t_IVCAC_N_0_5(minorizero(l)) = [];
     end 
     clear minorizero

      % Tempo IVC-min before AC
     for l = length(t_IVCminAC_G_0_5):-1:1
         if isnan(t_IVCminAC_G_0_5(l))
             t_IVCminAC_G_0_5(l) = [];
         end
     end
     minorizero = find(t_IVCminAC_G_0_5<0);
     for l = length(minorizero):-1:1
         t_IVCminAC_G_0_5(minorizero(l)) = [];
     end 
     clear minorizero

     for l = length(t_IVCminAC_N_0_5):-1:1
         if isnan(t_IVCminAC_N_0_5(l))
             t_IVCminAC_N_0_5(l) = [];
         end
     end
     minorizero = find(t_IVCminAC_N_0_5<0);
     for l = length(minorizero):-1:1
         t_IVCminAC_N_0_5(minorizero(l)) = [];
     end 
     clear minorizero
    
     % T IVC-MC
     for l = length(t_IVCMC_G0):-1:1
         if isnan(t_IVCMC_G0(l))
             t_IVCMC_G0(l) = [];
         end
     end
     minorizero = find(t_IVCMC_G0<0);
     for l = length(minorizero):-1:1
         t_IVCMC_G0(minorizero(l)) = [];
     end 
     clear minorizero

     for l = length(t_IVCMC_N0):-1:1
         if isnan(t_IVCMC_N0(l))
             t_IVCMC_N0(l) = [];
         end
     end
     minorizero = find(t_IVCMC_N0<0);
     for l = length(minorizero):-1:1
         t_IVCMC_N0(minorizero(l)) = [];
     end 
     clear minorizero

     for l = length(t_IVCMC_G_0_5):-1:1
         if isnan(t_IVCMC_G_0_5(l))
             t_IVCMC_G_0_5(l) = [];
         end
     end
     minorizero = find(t_IVCMC_G_0_5<0);
     for l = length(minorizero):-1:1
         t_IVCMC_G_0_5(minorizero(l)) = [];
     end 
     clear minorizero

     for l = length(t_IVCMC_N_0_5):-1:1
         if isnan(t_IVCMC_N_0_5(l))
             t_IVCMC_N_0_5(l) = [];
         end
     end
     minorizero = find(t_IVCMC_N_0_5<0);
     for l = length(minorizero):-1:1
         t_IVCMC_N_0_5(minorizero(l)) = [];
     end 

    % PEP
    for l = length(PEP_G_0_5):-1:1
         if isnan(PEP_G_0_5(l))
             PEP_G_0_5(l) = [];
         end
     end
     minorizero = find(PEP_G_0_5<0);
     for l = length(minorizero):-1:1
         PEP_G_0_5(minorizero(l)) = [];
     end 
     clear minorizero

     for l = length(PEP_N_0_5):-1:1
         if isnan(PEP_N_0_5(l))
             PEP_N_0_5(l) = [];
         end
     end
     minorizero = find(PEP_N_0_5<0);
     for l = length(minorizero):-1:1
         PEP_N_0_5(minorizero(l)) = [];
     end 
    % R_AO
     for l = length(R_AO_G_0_5):-1:1
         if isnan(R_AO_G_0_5(l))
             R_AO_G_0_5(l) = [];
         end
     end
     minorizero = find(R_AO_G_0_5<0);
     for l = length(minorizero):-1:1
         R_AO_G_0_5(minorizero(l)) = [];
     end 
     clear minorizero

     for l = length(R_AO_N_0_5):-1:1
         if isnan(R_AO_N_0_5(l))
             R_AO_N_0_5(l) = [];
         end
     end
     minorizero = find(R_AO_N_0_5<0);
     for l = length(minorizero):-1:1
         R_AO_N_0_5(minorizero(l)) = [];
     end 
% R_AC1
     for l = length(R_AC1_G_0_5):-1:1
         if isnan(R_AC1_G_0_5(l))
             R_AC1_G_0_5(l) = [];
         end
     end
     minorizero = find(R_AC1_G_0_5<0);
     for l = length(minorizero):-1:1
         R_AC1_G_0_5(minorizero(l)) = [];
     end 
     clear minorizero

     for l = length(R_AC1_N_0_5):-1:1
         if isnan(R_AC1_N_0_5(l))
             R_AC1_N_0_5(l) = [];
         end
     end
     minorizero = find(R_AC1_N_0_5<0);
     for l = length(minorizero):-1:1
         R_AC1_N_0_5(minorizero(l)) = [];
     end 

     % R_MC1
     for l = length(R_MC1_G_0_5):-1:1
         if isnan(R_MC1_G_0_5(l))
             R_MC1_G_0_5(l) = [];
         end
     end
     minorizero = find(R_MC1_G_0_5<0);
     for l = length(minorizero):-1:1
         R_MC1_G_0_5(minorizero(l)) = [];
     end 
     clear minorizero

     for l = length(R_MC1_N_0_5):-1:1
         if isnan(R_MC1_N_0_5(l))
             R_MC1_N_0_5(l) = [];
         end
     end
     minorizero = find(R_MC1_N_0_5<0);
     for l = length(minorizero):-1:1
         R_MC1_N_0_5(minorizero(l)) = [];
     end 
     
    %% Calcolo quanti tag tag 0 e tag 5 ho consecutivi (lunghezza massima e mediana di ogni finestra)
    % TUTTA LA GIORNATA
% Durante tutta la giornata quanti tag 0 o 5 ho    
[median_0,max_0,M_0] = maxmedian(0,R);
[median_5,max_5,M_5] = maxmedian(5,R);
% Durante tutta la giornata 0 e 5 INSIEME
 R_pp = R;
for i = 1:length(R_pp)
    if R_pp(i,3) == 5
        R_pp(i,3) = 0;
    end
end 
[median_05,max_05,M_05] = maxmedian(0,R_pp);
% Durante la giornata 0,4 e 5 INSIEME
R_pp1 = R_pp;
for i = 1:length(R_pp1)
    if R_pp1(i,3) == 4
        R_pp1(i,3) = 0;
    end
end 
[median_054,max_054,M_054] = maxmedian(0,R_pp1);

    % GIORNO
% Durante il giorno quanti tag 0 o 5 ho 
[median_0_G,max_0_G,M_0_G] = maxmedian(0,R_giorno);
[median_5_G,max_5_G,M_5_G] = maxmedian(5,R_giorno);
% Durante il giorno 0 e 5 INSIEME
R_giorno05 = R_pp(Giorno,:);
[median_05_G,max_05_G,M_05_G] = maxmedian(0,R_giorno05);
% Durante il giorno 0,4 e 5 INSIEME
R_giorno054 = R_pp1(Giorno,:);
[median_054_G,max_054_G,M_054_G] = maxmedian(0,R_giorno054);
    % NOTTE
 % Durante la notte quanti tag 0 o 5 ho    
[median_0_N,max_0_N,M_0_N] = maxmedian(0,R_notte);
[median_5_N,max_5_N,M_5_N] = maxmedian(5,R_notte);
% Durante la notte 0 e 5 INSIEME
R_notte05 = R_pp(Notte,:);
[median_05_N,max_05_N,M_05_N] = maxmedian(0,R_notte05);
% Durante la notte 0,4 e 5 INSIEME
R_notte054 = R_pp1(Notte,:);
[median_054_N,max_054_N,M_054_N] = maxmedian(0,R_notte054);


%%
    figure()
    subplot(131), histogram(M_0(:,2)),title('Number of beats - Tag 0'),xlabel('[Num of consecutive heartbeats]'),ylabel('[Occurrence]')
    subplot(132), histogram(M_5(:,2)),title('Number of beats - Tag 5'),xlabel('[Num of consecutive heartbeats]'),ylabel('[Occurrence]')
    subplot(133), histogram(M_05(:,2)),title('Number of beats - Tag 0 & 5'),xlabel('[Num of consecutive heartbeats]'),ylabel('[Occurrence]')
    sgtitle('Histogram All day')
    median_0 
    median_5
    median_05
    max_0
    max_5
    max_05
    pause
    figure()
    subplot(131), histogram(M_0_G(:,2)),title('Number of beats - Tag 0'),xlabel('[Num of consecutive heartbeats]'),ylabel('[Occurrence]')
    subplot(132), histogram(M_5_G(:,2)),title('Number of beats - Tag 5'),xlabel('[Num of consecutive heartbeats]'),ylabel('[Occurrence]')
    subplot(133), histogram(M_05_G(:,2)),title('Number of beats - Tag 0 & 5'),xlabel('[Num of consecutive heartbeats]'),ylabel('[Occurrence]')
    sgtitle('Histogram Day')
    median_0_G 
    median_5_G
    median_05_G
    max_0_G
    max_5_G
    max_05_G
    pause
    figure()
    subplot(131), histogram(M_0_N(:,2)),title('Number of beats - Tag 0'),xlabel('[Num of consecutive heartbeats]'),ylabel('[Occurrence]')
    subplot(132), histogram(M_5_N(:,2)),title('Number of beats - Tag 5'),xlabel('[Num of consecutive heartbeats]'),ylabel('[Occurrence]')
    subplot(133), histogram(M_05_N(:,2)),title('Number of beats - Tag 0 & 5'),xlabel('[Num of consecutive heartbeats]'),ylabel('[Occurrence]')
    sgtitle('Histogram Night')
    median_0_N 
    median_5_N
    median_05_N
    max_0_N
    max_5_N
    max_05_N
    pause
%%
     % Histogram - Amplitudes - TAG 0
     figure()
     subplot(321),histogram(amp_IVCAO_G0); hold on; histogram(amp_IVCAO_N0),xlabel('[mV]'),title('Amplitude IVC-AO'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(322),histogram(amp_IVCAC_G0); hold on; histogram(amp_IVCAC_N0),xlabel('[mV]'),title('Amplitude IVC-AC'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(323),histogram(amp_IVCMC_G0); hold on; histogram(amp_IVCMC_N0),xlabel('[mV]'),title('Amplitude IVC-MC'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(324),histogram(amp_IVCRE_G0); hold on; histogram(amp_IVCRE_N0),xlabel('[mV]'),title('Amplitude IVC-RE'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(325),histogram(amp_IVCminAC_G0); hold on; histogram(amp_IVCminAC_N0),xlabel('[mV]'),title('Amplitude IVC-minimum before AC'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(326),histogram(amp_IVCminAORE_G0); hold on; histogram(amp_IVCminAORE_N0),xlabel('[mV]'),title('Amplitude IVC-minimum between AO and RE'),legend('Day','Night'),ylabel('[Occurrence]')
     sgtitle('Amplitudes - Tag 0')

    % Histogram - Time intervals - TAG 0
     figure()
     subplot(321),histogram(t_IVCAO_G0.*1000); hold on; histogram(t_IVCAO_N0.*1000),xlabel('[ms]'),title('Time interval IVC-AO'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(322),histogram(t_IVCAC_G0.*1000); hold on; histogram(t_IVCAC_N0.*1000),xlabel('[ms]'),title('Time interval IVC-AC'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(323),histogram(t_IVCMC_G0.*1000); hold on; histogram(t_IVCMC_N0.*1000),xlabel('[ms]'),title('Time interval IVC-MC'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(324),histogram(t_IVCRE_G0.*1000); hold on; histogram(t_IVCRE_N0.*1000),xlabel('[ms]'),title('Time interval IVC-RE'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(325),histogram(t_IVCminAC_G0.*1000); hold on; histogram(t_IVCminAC_N0.*1000),xlabel('[ms]'),title('Time interval IVC-minimum before AC'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(326),histogram(t_IVCminAORE_G0.*1000); hold on; histogram(t_IVCminAORE_N0.*1000),xlabel('[ms]'),title('Time interval IVC-minimum between AO and RE'),legend('Day','Night'),ylabel('[Occurrence]')
     sgtitle('Time Intervals - Tag 0')
 
     % Histogram - Slopes - TAG 0
     figure()
     subplot(131),histogram(slope_IVCAO_G0); hold on; histogram(slope_IVCAO_N0),title('Slope IVC-AO'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(132),histogram(slope_minAORERE_G0); hold on; histogram(slope_minAORERE_N0),title('Slope minimum between AO and RE-RE'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(133),histogram(slope_minACAC_G0); hold on; histogram(slope_minACAC_N0),title('Slope minimum before AC-AC'),legend('Day','Night'),ylabel('[Occurrence]')
     sgtitle('Slopes - Tag 0')
 
     % Histogram - Other parameters1 - TAG 0
     figure()
     subplot(221),histogram(LVET_G0.*1000); hold on; histogram(LVET_N0.*1000),xlabel('[ms]'),title('Left Ventricular Ejection Time'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(222),histogram(QS2_G0.*1000); hold on; histogram(QS2_N0.*1000),xlabel('[ms]'),title('QS2'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(223),histogram(QT_G0.*1000); hold on; histogram(QT_N0.*1000),xlabel('[ms]'),title('QT'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(224),histogram(QTc_G0.*1000); hold on; histogram(QTc_N0.*1000),xlabel('[ms]'),title('Corrected QT'),legend('Day','Night'),ylabel('[Occurrence]')
     sgtitle('Other Parameters - Tag 0')

      % Histogram - Other parameters2 - TAG 0
     figure()
     subplot(221),histogram(PEP_G0.*1000); hold on; histogram(PEP_N0.*1000),xlabel('[ms]'),title('Pre Ejection Period'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(222),histogram(R_AO_G0.*1000); hold on; histogram(R_AO_N0.*1000),xlabel('[ms]'),title('Time Interval R-AO'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(223),histogram(R_AC1_G0.*1000); hold on; histogram(R_AC1_N0.*1000),xlabel('[ms]'),title('Time Interval R-AC'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(224),histogram(R_MC1_G0.*1000); hold on; histogram(R_MC1_N0.*1000),xlabel('[ms]'),title('Time Interval R-MC'),legend('Day','Night'),ylabel('[Occurrence]')
     sgtitle('Other Parameters - Tag 0')

     
      %% Histogram - Amplitudes - TAG 0 E 5
     figure()
     subplot(321),histogram(amp_IVCAO_G_0_5); hold on; histogram(amp_IVCAO_N_0_5),xlabel('[mV]'),title('Amplitude IVC-AO'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(322),histogram(amp_IVCAC_G_0_5); hold on; histogram(amp_IVCAC_N_0_5),xlabel('[mV]'),title('Amplitude IVC-AC'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(323),histogram(amp_IVCMC_G_0_5); hold on; histogram(amp_IVCMC_N_0_5),xlabel('[mV]'),title('Amplitude IVC-MC'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(324),histogram(amp_IVCRE_G_0_5); hold on; histogram(amp_IVCRE_N_0_5),xlabel('[mV]'),title('Amplitude IVC-RE'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(325),histogram(amp_IVCminAC_G_0_5); hold on; histogram(amp_IVCminAC_N_0_5),xlabel('[mV]'),title('Amplitude IVC-minimum before AC'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(326),histogram(amp_IVCminAORE_G_0_5); hold on; histogram(amp_IVCminAORE_N_0_5),xlabel('[mV]'),title('Amplitude IVC-minimum between AO and RE'),legend('Day','Night'),ylabel('[Occurrence]')
     sgtitle('Amplitudes - Tag 0 & 5')

     % Histogram - Time intervals  - TAG 0 E 5 
     figure()
     subplot(321),histogram(t_IVCAO_G_0_5.*1000); hold on; histogram(t_IVCAO_N_0_5.*1000),xlabel(('[ms]')),title('Time interval IVC-AO'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(322),histogram(t_IVCAC_G_0_5.*1000); hold on; histogram(t_IVCAC_N_0_5.*1000),xlabel(('[ms]')),title('Time interval IVC-AC'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(323),histogram(t_IVCMC_G_0_5.*1000); hold on; histogram(t_IVCMC_N_0_5.*1000),xlabel(('[ms]')),title('Time interval IVC-MC'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(324),histogram(t_IVCRE_G_0_5.*1000); hold on; histogram(t_IVCRE_N_0_5.*1000),xlabel(('[ms]')),title('Time interval IVC-RE'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(325),histogram(t_IVCminAC_G_0_5.*1000); hold on; histogram(t_IVCminAC_N_0_5.*1000),xlabel(('[ms]')),title('Time interval IVC-minimum before AC'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(326),histogram(t_IVCminAORE_G_0_5.*1000); hold on; histogram(t_IVCminAORE_N_0_5.*1000),xlabel(('[ms]')),title('Time interval IVC-minimum between AO and RE'),legend('Day','Night'),ylabel('[Occurrence]')
     sgtitle('Time Intervals - Tag 0 & 5')
 
     % Histogram - Slope - TAG 0 E 5
     figure()
     subplot(131),histogram(slope_IVCAO_G_0_5); hold on; histogram(slope_IVCAO_N_0_5),title('Slope IVC-AO'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(132),histogram(slope_minAORERE_G_0_5); hold on; histogram(slope_minAORERE_N_0_5),title('Slope minimum between AO and RE-RE'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(133),histogram(slope_minACAC_G_0_5); hold on; histogram(slope_minACAC_N_0_5),title('Slope minimum before AC-AC'),legend('Day','Night'),ylabel('[Occurrence]')
     sgtitle('Slopes - Tag 0 & 5')
 
     % Histogram - Other parameters1 - TAG 0 E 5
     figure()
     subplot(221),histogram(LVET_G_0_5.*1000); hold on; histogram(LVET_N_0_5.*1000),xlabel('[ms]'),title('Left Ventricular Ejection Time'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(222),histogram(QS2_G_0_5.*1000); hold on; histogram(QS2_N_0_5.*1000),xlabel('[ms]'),title('QS2'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(223),histogram(QT_G_0_5.*1000); hold on; histogram(QT_N_0_5.*1000),xlabel('[ms]'),title('QT'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(224),histogram(QTc_G_0_5.*1000); hold on; histogram(QTc_N_0_5.*1000),xlabel('[ms]'),title('Corrected QT'),legend('Day','Night'),ylabel('[Occurrence]')
     sgtitle('Other Parameters - Tag 0 & 5')

       % Histogram - Other parameters2 - TAG 0 E 5
     figure()
     subplot(221),histogram(PEP_G_0_5.*1000); hold on; histogram(PEP_N_0_5.*1000),xlabel('[ms]'),title('Pre Ejection Period'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(222),histogram(R_AO_G_0_5.*1000); hold on; histogram(R_AO_N_0_5.*1000),xlabel('[ms]'),title('Time Interval R-AO'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(223),histogram(R_AC1_G_0_5.*1000); hold on; histogram(R_AC1_N_0_5.*1000),xlabel('[ms]'),title('Time Interval R-AC'),legend('Day','Night'),ylabel('[Occurrence]')
     subplot(224),histogram(R_MC1_G_0_5.*1000); hold on; histogram(R_MC1_N_0_5.*1000),xlabel('[ms]'),title('Time Interval R-MC'),legend('Day','Night'),ylabel('[Occurrence]')
     sgtitle('Other Parameters - Tag 0 & 5')

 end 
%%

     % Histogram - Amplitudes - TAG 0
     figure()
     subplot(311),histogram(amp_IVCAO_G0); hold on; histogram(amp_IVCAO_N0),xlabel(['mV']),title('Amplitude IVC-AO'),legend('Day','Night')
     subplot(312),histogram(amp_IVCAC_G0); hold on; histogram(amp_IVCAC_N0),xlabel(['mV']),title('Amplitude IVC-AC'),legend('Day','Night')
     subplot(313),histogram(amp_IVCMC_G0); hold on; histogram(amp_IVCMC_N0),xlabel(['mV']),title('Amplitude IVC-MC'),legend('Day','Night')
     sgtitle('Amplitudes - Tag 0')
     figure()
     subplot(311),histogram(amp_IVCRE_G0); hold on; histogram(amp_IVCRE_N0),xlabel(['mV']),title('Amplitude IVC-RE'),legend('Day','Night')
     subplot(312),histogram(amp_IVCminAC_G0); hold on; histogram(amp_IVCminAC_N0),xlabel(['mV']),title('Amplitude IVC-minimum before AC'),legend('Day','Night')
     subplot(313),histogram(amp_IVCminAORE_G0); hold on; histogram(amp_IVCminAORE_N0),xlabel(['mV']),title('Amplitude IVC-minimum between AO and RE'),legend('Day','Night')
     sgtitle('Amplitudes - Tag 0')

%% 
% Histogram - Time intervals - TAG 0
     figure()
     subplot(311),histogram(t_IVCAO_G0.*1000); hold on; histogram(t_IVCAO_N0.*1000),xlabel(['ms']),title('Time interval IVC-AO'),legend('Day','Night')
     subplot(312),histogram(t_IVCAC_G0.*1000); hold on; histogram(t_IVCAC_N0.*1000),xlabel(['ms']),title('Time interval IVC-AC'),legend('Day','Night')
     subplot(313),histogram(t_IVCMC_G0.*1000); hold on; histogram(t_IVCMC_N0.*1000),xlabel(['ms']),title('Time interval IVC-MC'),legend('Day','Night')
     sgtitle('Time Intervals - Tag 0')
     figure()
     subplot(311),histogram(t_IVCRE_G0.*1000); hold on; histogram(t_IVCRE_N0.*1000),xlabel(['ms']),title('Time interval IVC-RE'),legend('Day','Night')
     subplot(312),histogram(t_IVCminAC_G0.*1000); hold on; histogram(t_IVCminAC_N0.*1000),xlabel(['ms']),title('Time interval IVC-minimum before AC'),legend('Day','Night')
     subplot(313),histogram(t_IVCminAORE_G0.*1000); hold on; histogram(t_IVCminAORE_N0.*1000),xlabel(['ms']),title('Time interval IVC-minimum between AO and RE'),legend('Day','Night')
     sgtitle('Time Intervals - Tag 0')
 %%
    figure()
     subplot(311),histogram(amp_IVCAO_G_0_5); hold on; histogram(amp_IVCAO_N_0_5),xlabel(['mV']),title('Amplitude IVC-AO'),legend('Day','Night')
     subplot(312),histogram(amp_IVCAC_G_0_5); hold on; histogram(amp_IVCAC_N_0_5),xlabel(['mV']),title('Amplitude IVC-AC'),legend('Day','Night')
     subplot(313),histogram(amp_IVCMC_G_0_5); hold on; histogram(amp_IVCMC_N_0_5),xlabel(['mV']),title('Amplitude IVC-MC'),legend('Day','Night')
     sgtitle('Amplitudes - Tag 0 & 5')
     figure()
     subplot(311),histogram(amp_IVCRE_G_0_5); hold on; histogram(amp_IVCRE_N_0_5),xlabel(['mV']),title('Amplitude IVC-RE'),legend('Day','Night')
     subplot(312),histogram(amp_IVCminAC_G_0_5); hold on; histogram(amp_IVCminAC_N_0_5),xlabel(['mV']),title('Amplitude IVC-minimum before AC'),legend('Day','Night')
     subplot(313),histogram(amp_IVCminAORE_G_0_5); hold on; histogram(amp_IVCminAORE_N_0_5),xlabel(['mV']),title('Amplitude IVC-minimum between AO and RE'),legend('Day','Night')
     sgtitle('Amplitudes - Tag 0 & 5')

     %%

     % Histogram - Time intervals  - TAG 0 E 5 
     figure()
     subplot(311),histogram(t_IVCAO_G_0_5.*1000); hold on; histogram(t_IVCAO_N_0_5.*1000),xlabel(['ms']),title('Time interval IVC-AO'),legend('Day','Night')
     subplot(312),histogram(t_IVCAC_G_0_5.*1000); hold on; histogram(t_IVCAC_N_0_5.*1000),xlabel(['ms']),title('Time interval IVC-AC'),legend('Day','Night')
     subplot(313),histogram(t_IVCMC_G_0_5.*1000); hold on; histogram(t_IVCMC_N_0_5.*1000),xlabel(['ms']),title('Time interval IVC-MC'),legend('Day','Night')
     sgtitle('Time Intervals - Tag 0 & 5')
     figure()
     subplot(311),histogram(t_IVCRE_G_0_5.*1000); hold on; histogram(t_IVCRE_N_0_5.*1000),xlabel(['ms']),title('Time interval IVC-RE'),legend('Day','Night')
     subplot(312),histogram(t_IVCminAC_G_0_5.*1000); hold on; histogram(t_IVCminAC_N_0_5.*1000),xlabel(['ms']),title('Time interval IVC-minimum before AC'),legend('Day','Night')
     subplot(313),histogram(t_IVCminAORE_G_0_5.*1000); hold on; histogram(t_IVCminAORE_N_0_5.*1000),xlabel(['ms']),title('Time interval IVC-minimum between AO and RE'),legend('Day','Night')
     sgtitle('Time Intervals - Tag 0 & 5')

    %%
%     R_MC1(:,2)=[];
%     R_MCsec(:,2)=[];
    R_MC1(:,2) = R(:,3);
    R_MCsec(:,2) = R(:,3);
    R_MC1_contr = R_MC1;
    R_MCsec_contr = R_MCsec;
    for i = length(R):-1:1
        if (R(i,3) == 1 || R(i,3) == 2 ||R(i,3) == 3 || R(i,3) == 6 || R(i,3) == 7 || R(i,3) == 8)
            R_MC1(i,:) = [];
            R_MCsec(i,:)=[];
        end
    end 
    figure()
    histogram(R_MC1(:,1))
    figure()
    histogram(R_MCsec(:,1))

    magg0 = find(R_MC1(:,1)>=0);
    length(magg0)
    POSSORECUPERARE = find(R_MCsec(:,1)>=-0.04);
    length(POSSORECUPERARE)