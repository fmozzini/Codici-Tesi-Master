%% Extract features from SCG
function [tIVCAO,tIVCAC,ampIVCAO,ampIVCAC,slopeIVCAO,LVET,QS2,QT,QTc,tIVCMC,tIVCRE,tIVCminAORE,tIVCminAC,ampIVCMC,ampIVCRE,ampIVCminAORE,ampIVCminAC,...
    slopeminAORERE,slopeminACAC,RdivT] = extractfeatures(fs_scg,fs_ecg,posAO,ampAO,posIVC,ampIVC,posAC,ampAC,RR,posQ_SCG,posfineT,ampT,posQ_ECG,posMC,...
    ampMC,posRE,ampRE,posminAORE,ampminAORE,posminAC,ampminAC,ampR)
    % Temporal distances
    tIVCAO = (posAO-posIVC)/fs_scg;
    tIVCAC = (posAC-posIVC)/fs_scg;
    tIVCMC = (posIVC-posMC)/fs_scg;
    tIVCRE = (posRE-posIVC)/fs_scg;
    tIVCminAORE = (posminAORE-posIVC)/fs_scg;
    tIVCminAC = (posminAC - posIVC)/fs_scg;
    % Amplitudes
    ampIVCAO = ampAO-ampIVC;
    ampIVCAC = ampAC-ampIVC;
    ampIVCMC = ampMC-ampIVC;
    ampIVCRE = ampRE-ampIVC;
    ampIVCminAORE = ampminAORE-ampIVC;
    ampIVCminAC = ampminAC-ampIVC;
    % Slopes
    slopeIVCAO = ampIVCAO/tIVCAO;
    slopeminAORERE = (ampRE-ampminAORE)/(posRE-posminAORE);
    slopeminACAC = (ampAC-ampminAC)/(posAC-posminAC);
       %% LVET, QS2, QT e QTc
    LVET = (posAC-posAO)/fs_scg; % [s]
    QS2 = (posAC-posQ_SCG)/fs_scg; % [s]
    QT = (posfineT-posQ_ECG)./fs_ecg;
    QTc = (QT/(sqrt(RR)))./fs_ecg; % [s]
    % Altri parametri
    RdivT = ampR/ampT;
   
end 