%% Codice per fare Cosinor Analysis e salvare direttamente tutti i dati in un file csv
clear all
close all
clc
%% 
folderPAR = 'C:\Users\feder\Desktop\Tesi\Data\Parameters SCG_10SEC';
listPAR = dir(folderPAR);
listPAR(1) = [];
listPAR(1) = [];

 addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Parameters SCG_10SEC'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Codes\
 %%
 for m = 1:1
    FOLDERPAR = fullfile(listPAR(m).folder, listPAR(m).name)
    file = dir(FOLDERPAR);
    name = file.name;
    load(name)

    Amplitudes = [amp_IVCAO(:,1) amp_IVCMC(:,1) amp_IVCRE(:,1) amp_IVCAC(:,1) amp_IVCminAORE(:,1) amp_IVCminAC(:,1)];
    Amplitudes_name = ["amp_IVCAO" "amp_IVCMC" "amp_IVCRE" "amp_IVCAC" "amp_IVCminAORE" "amp_IVCminAC"];
    Temporals = [t_IVCAO(:,1) t_IVCMC(:,1) t_IVCRE(:,1) t_IVCAC(:,1) t_IVCminAORE(:,1) t_IVCminAC(:,1)...
        LVET(:,1) QS2(:,1) QT(:,1) QTc(:,1) PEP(:,1) R_AO(:,1) R_AC1(:,1) R_MC1(:,1)];
    Temporals_name = ["t_IVCAO" "t_IVCMC" "t_IVCRE" "t_IVCAC" "t_IVCminAORE" "t_IVCminAC"...
        "LVET" "QS2" "QT" "QTc" "PEP" "R_AO" "R_AC1" "R_MC1"];
    Slopes = [slope_IVCAO(:,1) slope_minACAC(:,1) slope_minAORERE(:,1)];
    Slopes_name = ["slope_IVCAO" "slope_minACAC" "slope_minAORERE"];
    Parameters = [Amplitudes Temporals Slopes];
    Parameters_name = [Amplitudes_name Temporals_name Slopes_name];
    Results = zeros(size(Parameters,2),6);
    for c = 1:length(Parameters)
        try
        [MESOR, OA, phi, pval] = run_cosinor_fittedcurve(Parameters(:,c),Parameters_name(c));
        pause
        close all
        catch ME 
            continue; 
            % se c'è un errore va al parametro successivo, per quella riga
            % avrò poi tutti 0 (a parte il nome del parametro)
        end 
       
        % Riesco a trovare un modo per tenere conto di quando mi da
        % l'errore sulla confidence?
        Results(c,:) = [MESOR OA phi pval MESOR/OA OA/MESOR];
        Results_name = [char(Parameters_name') num2str(Results)];
    end
    name = erase(name,"Parameters SCG-ECG_FILT-")
% % %     save(['C:\Users\feder\Desktop\Tesi\Data\PostProc PT 1\' 'PostProc PT 1-' name], 'QRS','HR_min','HR_5min')
    save(['C:\Users\feder\Desktop\Tesi\Data\Cosinor Analysis\' 'Cosinor-' name], 'Results','Results_name')
    
 end 



