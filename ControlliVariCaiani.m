%% Conto quanti picchi per tipologia
n_picchi0 = 0;
n_picchi1 = 0;
n_picchi2 = 0;
n_picchipiudi3 = 0;
for i = 1:length(n_picchi)
    if n_picchi(i) == 0
        n_picchi0 = n_picchi0 + 1;
    elseif  n_picchi(i) == 1
        n_picchi1 = n_picchi1 + 1;
    elseif n_picchi(i) == 2
        n_picchi2 = n_picchi2 + 1;
    else 
        n_picchipiudi3 = n_picchipiudi3 + 1;        
    end
end 
%% Dei picchi 0, quanti sono durante il giorno e quanti durante la notte?
countnotte = 0;
countgiorno = 0;
for i = 1:length(R)-1
    if n_picchi(i) == 0 && GN(i) == 0
        countnotte = countnotte+1;
    elseif n_picchi(i) == 0 && GN(i) == 1
        countgiorno = countgiorno+1;
    end 
end 
%% Plotto per capire un po' cosa succede
fs_ECG = 1024;
fs_SCG = 64;
window_SCG = 0.2*fs_SCG;
window_ECG = 0.2*fs_ECG;
count = 0;
for i = 1:length(R)
    if n_picchi(i) == 0 && GN(i) == 1
        count = count+1;
        qrs1 = (R(i,1)/1024)*64;
        qrs2 = (R(i+1,1)/1024)*64;
        finestrabattito_ECGF = ECG_filt(R(i,1)-window_ECG:R(i+1,1)-window_ECG)';
        finestrabattito_SCG = Acc_z(qrs1-window_SCG:qrs2-window_SCG)';
        figure()
        set(gcf, 'WindowState', 'maximized');
        a = subplot(211); plot((R(i,1)-window_ECG:R(i+1,1)-window_ECG)./1024,finestrabattito_ECGF),xlabel('[s]'),title('Ecg')
        b = subplot(212); plot((qrs1-window_SCG:qrs2-window_SCG)./64,finestrabattito_SCG),xlabel('[s]'),title('SCG')
        sgtitle(i)
        pause
        close all
    end 
end 
%% Conto il numero di righe consecutive per tag0 e tag 5 (quanti battiti consecutivi ho durante il giorno e durante la notte)
% Tag 0, Notte
Maxduratanotte_0 = 0;
duratanotte_0 = 0;
for i = 2:length(R)-1
    if R(i,3) == 0 && GN(i) == 0 %notte
        if R(i-1,3) == 0 && GN(i-1) == 0
            duratanotte_0 = duratanotte_0+1;
            if duratanotte_0 > Maxduratanotte_0
                Maxduratanotte_0 = duratanotte_0;
            end 
        else
            duratanotte_0 = 0;
        end 
    end
end 

% Tag 0, Giorno
Maxduratagiorno_0 = 0;
duratagiorno_0 = 0;
for i = 2:length(R)-1
    if R(i,3) == 0 && GN(i) == 1 %giorno
        if R(i-1,3) == 0 && GN(i-1) == 1
            duratagiorno_0 = duratagiorno_0+1;
            if duratagiorno_0 > Maxduratagiorno_0
                Maxduratagiorno_0 = duratagiorno_0;
            end 
        else
            duratagiorno_0 = 0;
        end 
    end
end 
%%
% Tag 0 e 5, Giorno
Maxduratagiorno_05 = 0;
duratagiorno_05 = 0;
for i = 2:length(R)-1
    if R(i,3) == 0 || R(i,3) == 5
        if GN(i) == 1 %giorno
            if R(i-1,3) == 0 || R(i-1,3) == 5
                if GN(i-1) == 1
                    duratagiorno_05 = duratagiorno_05+1;
                    if duratagiorno_05 > Maxduratagiorno_05
                        Maxduratagiorno_05 = duratagiorno_05;
                    end 
                else
                    duratagiorno_05 = 0;
                end 
            end 
        end 
    end
end 

% Tag 0 e 5, Notte
Maxduratanotte_05 = 0;
duratanotte_05 = 0;
for i = 2:length(R)-1
    if R(i,3) == 0 || R(i,3) == 5
        if GN(i) == 0 %notte
            if R(i-1,3) == 0 || R(i-1,3) == 5
                if GN(i-1) == 0
                    duratanotte_05 = duratanotte_05+1;
                    if duratanotte_05 > Maxduratanotte_05
                        Maxduratanotte_05 = duratanotte_05;
                    end 
                else
                    duratanotte_05 = 0;
                end 
            end 
        end 
    end
end 
