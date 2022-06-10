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
    
    %% Calcolo quanti tag tag 0 e tag 5 ho consecutivi (lunghezza massima e mediana di ogni finestra)

    %% Controllo in generale (giorno e notte)
    M = zeros(length(R), 2);
    c = 0;
    ln = 1;
    val = R(1,3);
    for i = 2:length(R)
        c = c + 1;
        cur_val = R(i,3);
        if cur_val ~= val
            M(ln, 1) = val;
            M(ln, 2) = c;
            c = 0;
            val = cur_val;
            ln = ln+1;
        end
    end
    M(ln, 1) = val;
    M(ln, 2) = c;
    M(ln:end,:) = [];
    % Giorno e notte - tag 0
    posM_0 = find(M(:,1) == 0); M_0 = M(posM_0,:);
    median_0 = median(M_0(:,2)); max_0 = max(M_0(:,2));
    % Giorno e notte - tag 5 
    posM_5 = find(M(:,1) == 5); M_5 = M(posM_5,:);
    median_5 = median(M_5(:,2)); max_5 = max(M_5(:,2));
    %% Giorno e notte - tag 0 & tag 5
    R_pp = R;
    for i = 1:length(R_pp)
        if R_pp(i,3) == 5
            R_pp(i,3) = 0;
        end
    end 
    M05 = zeros(length(R_pp), 2);
    c05 = 0;
    ln05 = 1;
    val05 = R_pp(1,3);
    for i = 2:length(R_pp)
        c05 = c05 + 1;
        cur_val05 = R_pp(i,3);
        if cur_val05 ~= val05
            M05(ln05, 1) = val05;
            M05(ln05, 2) = c05;
            c05 = 0;
            val05 = cur_val05;
            ln05 = ln05+1;
        end
    end
    M05(ln05, 1) = val05;
    M05(ln05, 2) = c05;
    M05(ln05:end,:) = [];
    % Intero (giorno e notte insieme) tag 0 e 5 
    posM_05 = find(M05(:,1) == 0); M_05 = M05(posM_05,:);
    median_05 = median(M_05(:,2)); max_05 = max(M_05(:,2));

    %% Giorno
    M_giorno = zeros(length(R_giorno), 2);
    c_giorno = 0;
    ln_giorno = 1;
    val_giorno = R_giorno(1,3);
    for i = 2:length(R_giorno)
        c_giorno = c_giorno + 1;
        cur_val_giorno = R_giorno(i,3);
        if cur_val_giorno ~= val_giorno
            M_giorno(ln_giorno, 1) = val_giorno;
            M_giorno(ln_giorno, 2) = c_giorno;
            c_giorno = 0;
            val_giorno = cur_val_giorno;
            ln_giorno = ln_giorno+1;
        end
    end
    M_giorno(ln_giorno, 1) = val_giorno;
    M_giorno(ln_giorno, 2) = c_giorno;
    M_giorno(ln_giorno:end,:) = [];
    % GIORNO - tag 0 
    posM_0_G = find(M_giorno(:,1) == 0); M_0_G = M_giorno(posM_0_G,:); 
    median_0_G = median(M_0_G(:,2)); max_0_G = max(M_0_G(:,2));
    % GIORNO - tag 5 
    posM_5_G = find(M_giorno(:,1) == 5); M_5_G = M_giorno(posM_5_G,:);
    median_5_G = median(M_5_G(:,2)); max_5_G = max(M_5_G(:,2));

    %% Giorno - tag 0 & tag 5
    R_giorno05 = R_pp(Giorno,:);
    M_giorno05 = zeros(length(R_giorno05), 2);
    c_giorno05 = 0;
    ln_giorno05 = 1;
    val_giorno05 = R_giorno05(1,3);
    for i = 2:length(R_giorno05)
        c_giorno05 = c_giorno05 + 1;
        cur_val_giorno05 = R_giorno05(i,3);
        if cur_val_giorno05 ~= val_giorno05
            M_giorno05(ln_giorno05, 1) = val_giorno05;
            M_giorno05(ln_giorno05, 2) = c_giorno05;
            c_giorno05 = 0;
            val_giorno05 = cur_val_giorno05;
            ln_giorno05 = ln_giorno05+1;
        end
    end
    M_giorno05(ln_giorno05, 1) = val_giorno05;
    M_giorno05(ln_giorno05, 2) = c_giorno05;
    M_giorno05(ln_giorno05:end,:) = [];
    % Giorno -  tag 0 e 5 
    posM_05_G = find(M_giorno05(:,1) == 0); M_05_G = M_giorno05(posM_05_G,:);
    median_05_G = median(M_05_G(:,2)); max_05_G = max(M_05_G(:,2));

    %% Notte
    M_notte = zeros(length(R_notte), 2);
    c_notte = 0;
    ln_notte = 1;
    val_notte = R_notte(1,3);
    for i = 2:length(R_notte)
        c_notte = c_notte + 1;
        cur_val_notte = R_notte(i,3);
        if cur_val_notte ~= val_notte
            M_notte(ln_notte, 1) = val_notte;
            M_notte(ln_notte, 2) = c_notte;
            c_notte = 0;
            val_notte = cur_val_notte;
            ln_notte = ln_notte+1;
        end
    end
    M_notte(ln_notte, 1) = val_notte;
    M_notte(ln_notte, 2) = c_notte;
    M_notte(ln_notte:end,:) = [];
    % NOTTE - Tag 0 
    posM_0_N = find(M_notte(:,1) == 0); M_0_N = M_notte(posM_0_N,:); 
    median_0_N = median(M_0_N(:,2)); max_0_N = max(M_0_N(:,2));
    % NOTTE - Tag 5
    posM_5_N = find(M_notte(:,1) == 5); M_5_N = M_notte(posM_5_N,:);
    median_5_N = median(M_5_N(:,2)); max_5_N = max(M_5_N(:,2));

    %% Notte - tag 0 & tag 5
    R_notte05 = R_pp(Notte,:);
    M_notte05 = zeros(length(R_notte05), 2);
    c_notte05 = 0;
    ln_notte05 = 1;
    val_notte05 = R_notte05(1,3);
    for i = 2:length(R_notte05)
        c_notte05 = c_notte05 + 1;
        cur_val_notte05 = R_notte05(i,3);
        if cur_val_notte05 ~= val_notte05
            M_notte05(ln_notte05, 1) = val_notte05;
            M_notte05(ln_notte05, 2) = c_notte05;
            c_notte05 = 0;
            val_notte05 = cur_val_notte05;
            ln_notte05 = ln_notte05+1;
        end
    end
    M_notte05(ln_notte05, 1) = val_notte05;
    M_notte05(ln_notte05, 2) = c_notte05;
    M_notte05(ln_notte05:end,:) = [];
    % Night -  tag 0 e 5 
    posM_05_N = find(M_notte05(:,1) == 0); M_05_N = M_notte05(posM_05_N,:);
    median_05_N = median(M_05_N(:,2)); max_05_N = max(M_05_N(:,2));

    %% AGGIUNGO ANCHE IL TAG 4 (tag0,4,5)
     %% Giorno e notte - tag 0 & tag 5 % tag 5
    R_pp1 = R;
    for i = 1:length(R_pp1)
        if R_pp1(i,3) == 5
            R_pp1(i,3) = 0;
        elseif R_pp1(i,3) == 4
            R_pp1(i,3) = 0;
        end
    end 
    M054 = zeros(length(R_pp1), 2);
    c054 = 0;
    ln054 = 1;
    val054 = R_pp1(1,3);
    for i = 2:length(R_pp1)
        c054 = c054 + 1;
        cur_val054 = R_pp1(i,3);
        if cur_val054 ~= val054
            M054(ln054, 1) = val054;
            M054(ln054, 2) = c054;
            c054 = 0;
            val054 = cur_val054;
            ln054 = ln054+1;
        end
    end
    M054(ln054, 1) = val054;
    M054(ln054, 2) = c054;
    M054(ln054:end,:) = [];
    % Intero (giorno e notte insieme) tag 0 e 5 e 4
    posM_054 = find(M054(:,1) == 0); M_054 = M05(posM_054,:);
    median_054 = median(M_054(:,2)); max_054 = max(M_054(:,2));

    %% Giorno - tag 0 & tag 5 & tag 4
    R_giorno054 = R_pp1(Giorno,:);
    M_giorno054 = zeros(length(R_giorno054), 2);
    c_giorno054 = 0;
    ln_giorno054 = 1;
    val_giorno054 = R_giorno054(1,3);
    for i = 2:length(R_giorno054)
        c_giorno054 = c_giorno054 + 1;
        cur_val_giorno054 = R_giorno054(i,3);
        if cur_val_giorno054 ~= val_giorno054
            M_giorno054(ln_giorno054, 1) = val_giorno054;
            M_giorno054(ln_giorno054, 2) = c_giorno054;
            c_giorno054 = 0;
            val_giorno054 = cur_val_giorno054;
            ln_giorno054 = ln_giorno054+1;
        end
    end
    M_giorno054(ln_giorno054, 1) = val_giorno054;
    M_giorno054(ln_giorno054, 2) = c_giorno054;
    M_giorno054(ln_giorno054:end,:) = [];
    % Giorno -  tag 0 e 5 e 4
    posM_054_G = find(M_giorno054(:,1) == 0); M_054_G = M_giorno05(posM_054_G,:);
    median_054_G = median(M_054_G(:,2)); max_054_G = max(M_054_G(:,2));

    %% Notte - tag 0 & tag 5 & tag 4
    R_notte054 = R_pp1(Notte,:);
    M_notte054 = zeros(length(R_notte054), 2);
    c_notte054 = 0;
    ln_notte054 = 1;
    val_notte054 = R_notte054(1,3);
    for i = 2:length(R_notte054)
        c_notte054 = c_notte054 + 1;
        cur_val_notte054 = R_notte054(i,3);
        if cur_val_notte054 ~= val_notte054
            M_notte054(ln_notte054, 1) = val_notte054;
            M_notte054(ln_notte054, 2) = c_notte054;
            c_notte054 = 0;
            val_notte054 = cur_val_notte054;
            ln_notte054 = ln_notte054+1;
        end
    end
    M_notte054(ln_notte054, 1) = val_notte054;
    M_notte054(ln_notte054, 2) = c_notte054;
    M_notte054(ln_notte054:end,:) = [];
    % Night -  tag 0 e 5 e 4 
    posM_054_N = find(M_notte054(:,1) == 0); M_054_N = M_notte054(posM_054_N,:);
    median_054_N = median(M_054_N(:,2)); max_054_N = max(M_054_N(:,2));

%%
    figure()
    subplot(141), histogram(M_0(:,2)),title('Number of beats - Tag 0')
    subplot(142), histogram(M_5(:,2)),title('Number of beats - Tag 5')
    subplot(143), histogram(M_05(:,2)),title('Number of beats - Tag 0 & 5')
    subplot(144), histogram(M_054(:,2)),title('Number of beats - Tag 0 & 5 & 4')
    sgtitle('Histogram All day')
    median_0 
    median_5
    median_05
    median_054
    max_0
    max_5
    max_05
    max_054
    pause
    figure()
    subplot(141), histogram(M_0_G(:,2)),title('Number of beats - Tag 0')
    subplot(142), histogram(M_5_G(:,2)),title('Number of beats - Tag 5')
    subplot(143), histogram(M_05_G(:,2)),title('Number of beats - Tag 0 & 5')
    subplot(144), histogram(M_054_G(:,2)),title('Number of beats - Tag 0 & 5 & 4')
    sgtitle('Histogram Day')
    median_0_G 
    median_5_G
    median_05_G
    median_054_G
    max_0_G
    max_5_G
    max_05_G
    max_054_G
    pause
    figure()
    subplot(141), histogram(M_0_N(:,2)),title('Number of beats - Tag 0')
    subplot(142), histogram(M_5_N(:,2)),title('Number of beats - Tag 5')
    subplot(143), histogram(M_05_N(:,2)),title('Number of beats - Tag 0 & 5')
    subplot(144), histogram(M_054_N(:,2)),title('Number of beats - Tag 0 & 5 & 4')
    sgtitle('Histogram Night')
    median_0_N 
    median_5_N
    median_05_N
    median_054_N
    max_0_N
    max_5_N
    max_05_N
    max_054_N
    pause
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
     subplot(321),histogram(t_IVCAO_G0.*1000); hold on; histogram(t_IVCAO_N0.*1000),xlabel(['ms']),title('Time interval IVC-AO'),legend('Day','Night')
     subplot(322),histogram(t_IVCAC_G0.*1000); hold on; histogram(t_IVCAC_N0.*1000),xlabel(['ms']),title('Time interval IVC-AC'),legend('Day','Night')
     subplot(323),histogram(t_IVCMC_G0.*1000); hold on; histogram(t_IVCMC_N0.*1000),xlabel(['ms']),title('Time interval IVC-MC'),legend('Day','Night')
     subplot(324),histogram(t_IVCRE_G0.*1000); hold on; histogram(t_IVCRE_N0.*1000),xlabel(['ms']),title('Time interval IVC-RE'),legend('Day','Night')
     subplot(325),histogram(t_IVCminAC_G0.*1000); hold on; histogram(t_IVCminAC_N0.*1000),xlabel(['ms']),title('Time interval IVC-minimum before AC'),legend('Day','Night')
     subplot(326),histogram(t_IVCminAORE_G0.*1000); hold on; histogram(t_IVCminAORE_N0.*1000),xlabel(['ms']),title('Time interval IVC-minimum between AO and RE'),legend('Day','Night')
     sgtitle('Time Intervals - Tag 0')
 
     % Histogram - Slopes - TAG 0
     figure()
     subplot(131),histogram(slope_IVCAO_G0); hold on; histogram(slope_IVCAO_N0),title('Slope IVC-AO'),legend('Day','Night')
     subplot(132),histogram(slope_minAORERE_G0); hold on; histogram(slope_minAORERE_N0),title('Slope minimum between AO and RE-RE'),legend('Day','Night')
     subplot(133),histogram(slope_minACAC_G0); hold on; histogram(slope_minACAC_N0),title('Slope minimum before AC-AC'),legend('Day','Night')
     sgtitle('Slopes - Tag 0')
 
     % Histogram - Other parameters - TAG 0
     figure()
     subplot(221),histogram(LVET_G0.*1000); hold on; histogram(LVET_N0.*1000),xlabel(['ms']),title('Left Ventricular Ejection Time'),legend('Day','Night')
     subplot(222),histogram(QS2_G0.*1000); hold on; histogram(QS2_N0.*1000),xlabel(['ms']),title('QS2'),legend('Day','Night')
     subplot(223),histogram(QT_G0.*1000); hold on; histogram(QT_N0.*1000),xlabel(['ms']),title('QT'),legend('Day','Night')
     subplot(224),histogram(QTc_G0.*1000); hold on; histogram(QTc_N0.*1000),xlabel(['ms']),title('Corrected QT'),legend('Day','Night')
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
     subplot(321),histogram(t_IVCAO_G_0_5.*1000); hold on; histogram(t_IVCAO_N_0_5.*1000),xlabel(['ms']),title('Time interval IVC-AO'),legend('Day','Night')
     subplot(322),histogram(t_IVCAC_G_0_5.*1000); hold on; histogram(t_IVCAC_N_0_5.*1000),xlabel(['ms']),title('Time interval IVC-AC'),legend('Day','Night')
     subplot(323),histogram(t_IVCMC_G_0_5.*1000); hold on; histogram(t_IVCMC_N_0_5.*1000),xlabel(['ms']),title('Time interval IVC-MC'),legend('Day','Night')
     subplot(324),histogram(t_IVCRE_G_0_5.*1000); hold on; histogram(t_IVCRE_N_0_5.*1000),xlabel(['ms']),title('Time interval IVC-RE'),legend('Day','Night')
     subplot(325),histogram(t_IVCminAC_G_0_5.*1000); hold on; histogram(t_IVCminAC_N_0_5.*1000),xlabel(['ms']),title('Time interval IVC-minimum before AC'),legend('Day','Night')
     subplot(326),histogram(t_IVCminAORE_G_0_5.*1000); hold on; histogram(t_IVCminAORE_N_0_5.*1000),xlabel(['ms']),title('Time interval IVC-minimum between AO and RE'),legend('Day','Night')
     sgtitle('Time Intervals - Tag 0 & 5')
 
     % Histogram - Slope - TAG 0 E 5
     figure()
     subplot(131),histogram(slope_IVCAO_G_0_5); hold on; histogram(slope_IVCAO_N_0_5),title('Slope IVC-AO'),legend('Day','Night')
     subplot(132),histogram(slope_minAORERE_G_0_5); hold on; histogram(slope_minAORERE_N_0_5),title('Slope minimum between AO and RE-RE'),legend('Day','Night')
     subplot(133),histogram(slope_minACAC_G_0_5); hold on; histogram(slope_minACAC_N_0_5),title('Slope minimum before AC-AC'),legend('Day','Night')
     sgtitle('Slopes - Tag 0 & 5')
 
     % Histogram - Other parameters - TAG 0 E 5
     figure()
     subplot(221),histogram(LVET_G_0_5.*1000); hold on; histogram(LVET_N_0_5.*1000),xlabel(['ms']),title('Left Ventricular Ejection Time'),legend('Day','Night')
     subplot(222),histogram(QS2_G_0_5.*1000); hold on; histogram(QS2_N_0_5.*1000),xlabel(['ms']),title('QS2'),legend('Day','Night')
     subplot(223),histogram(QT_G_0_5.*1000); hold on; histogram(QT_N_0_5.*1000),xlabel(['ms']),title('QT'),legend('Day','Night')
     subplot(224),histogram(QTc_G_0_5.*1000); hold on; histogram(QTc_N_0_5.*1000),xlabel(['ms']),title('Corrected QT'),legend('Day','Night')
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

    
