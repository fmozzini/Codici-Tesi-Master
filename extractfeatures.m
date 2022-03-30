%% Extract features from SCG
function [tIVCAO,tIVCAC,ampIVCAO,ampIVCAC,slopeIVCAO] = extractfeatures(posAO,ampAO,posIVC,ampIVC,posAC,ampAC,fs)
    tIVCAO = (posAO-posIVC)/fs;
    tIVCAC = (posAC-posIVC)/fs;
    ampIVCAO = ampAO-ampIVC;
    ampIVCAC = ampAC-ampIVC;
    slopeIVCAO = ampIVCAO/tIVCAO;
end 