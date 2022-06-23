function [MESOR, OA, phi, pval] = run_cosinor_fittedcurve(vect)
% INPUT: 
% vect = series on which to perform the cosino analysis 

% OUTPUT:
% MESOR = midline value
% OA = day/night oscillation amplitude  
% phi = acrophase, temporal value at which the amplitude of the fitting sinusoid is maximal
% pval = p value of the Zero Amplitude Test: it is 1 if pval < 0.05. If 1 = presence of circadian rhythm 

%% Pre-process

% Interpolate nan
nanx = isnan(vect);
t    = 1:numel(vect); 
vect(nanx) = interp1(t(~nanx), vect(~nanx), t(nanx));

if isnan(vect(1,1))
    vect = fillmissing(vect,'next');    
end
if isnan(vect(end,1))
    vect = fillmissing(vect,'previous');    
end

%%

w=2*pi;
alpha=0.05;
t=[1:size(vect, 1 )]/size(vect, 1 );
t=t';
[MESOR, OA, phi, RSS, CI_Amp_min, CI_Amp_max, pval, ~]=cosinor(t,vect,w,alpha);

close
% pause
end

