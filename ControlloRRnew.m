%% Federica Mozzini - 946400
% Ricalcoliamo gli intervalli RR e guardiamo cosa succede in ECG quando RR
% è strano

clear all
close all
clc
%%
folderSCG = 'C:\Users\feder\Desktop\Tesi\Data\Picchi SCG - Acc z';
folderPT = 'C:\Users\feder\Desktop\Tesi\Data\Pan-Tompkins'; 
folderECG = 'C:\Users\feder\Desktop\Tesi\Data\Filtered ECG';

listPT = dir(folderPT);
listPT(1) = [];
listPT(1) = [];
listECG = dir(folderECG);
listECG(1) = [];
listECG(1) = [];
listECG(end) = [];
list = dir(folderSCG);
list(1) = [];
list(1) = [];
N = length(list);
% list(N-1) = [];
% list(N-1) = [];
% N = length(list)

 addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Filtered ECG'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Data\Pan-Tompkins\
 addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Picchi SCG - Acc z'\
 addpath 'C:\Users\feder\Desktop\Tesi'\Codes\

 %% 
 for m = 16:16
    FOLDERSCG = fullfile(list(m).folder, list(m).name)
    file = dir(FOLDERSCG);
    name = file.name;
    load(name)

    FOLDERPT = fullfile(listPT(m).folder, listPT(m).name)
    file = dir(FOLDERPT);
    name = file.name;
    load(name)

    FOLDERECG = fullfile(listECG(m).folder, listECG(m).name)
    file = dir(FOLDERECG);
    name = file.name;
    load(name)

    fs_ECG = 1024;
    fs_SCG = 64;
    window_SCG = 0.2*fs_SCG;
    window_ECG = 0.2*fs_ECG;
    R = zeros(peaksFORwindow(end-1),2);

    for i = 1:peaksFORwindow(end-1)
        R(i,:) = [qrs_I(i) qrs_AMP(i)];
    end
   
    RR = diff(R(:,1)); %indici, non s o ms (divido per 1024)
    RR_sec = RR./1024; % lo plotto in sec sulle y, 84163
    RR_tempo = R(1:end-1,1);
    RR_tempo = RR_tempo./1024;
%     plot(RR_tempo,RR_sec),title('RR interval PRE PROCESSING')

%% CERCHIAMO UN MODO PER TROVARE GLI OUTLIERS -> elimino i picchi se sono diversi dal 30% della media tra valore precedente e successivo
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6931967/ 
    eliminare = 0;
        for r = 1:length(RR_sec)-2 
            media = (RR_sec(r)+RR_sec(r+2))/2;
            diff = media - RR_sec(r+1);
            realdiff = diff/media;
            if (abs(realdiff)) >= 0.3
                eliminare = eliminare + 1;
                ELIMINARE(eliminare,:) = [RR_sec(r+1) RR_tempo(r+1) r+1];
            end 
        end 
    RR_pp = [RR_sec RR_tempo];
    R_post_processing = R;
    for n = length(ELIMINARE):-1:1
    %     R_post_processing(find(R_post_processing ==  ELIMINARE(n,2)),:) = [];
        R_post_processing(ELIMINARE(n,3),:) = [];
        RR_pp(ELIMINARE(n,3),:) = [];
    end 
   % ne ho eliminati 480 ma rimane qualcosina (pochi pochi outliers)
%    figure()
%    plot(RR_pp(:,2),RR_pp(:,1)),xlabel('[s]'),ylabel('RR[s]'),title('RR interval POST PROCESSING'); 


%     SALVO IL RISULTATO
    name = erase(name,"ECG_FILT-")
    save(['C:\Users\feder\Desktop\Tesi\Data\PostProc PT\' 'PostProc PT-' name],'R_post_processing','RR_pp','peaksFORwindow','HR_min','HR_5min')

 end 

        %% GRAFICO PER PRESENTAZIONE - TESI
figure()
plot(RR_tempo,RR_sec); hold on;plot(RR_pp(:,2),RR_pp(:,1),'r'),xlabel('[s]'),ylabel('RR[s]'),legend('Pre Processing','Post Processing'),...
    title('RR interval')

cinqueminuti_sec = 5*60;
RR = [RR_sec RR_tempo];
for r = length(RR):-1:1 
    if RR(r,2)<cinqueminuti_sec 
        RR(r,:) = [];
    end
end 
fine = RR(end,2)-cinqueminuti_sec;
conta = 0;
for r = length(RR):-1:1
    if RR(r,2)>fine
        conta = conta + 1;
        RR(r,:) = [];
    end 
end 

RR_pp1 = RR_pp;
for r = length(RR_pp1):-1:1 
    if RR_pp1(r,2)<cinqueminuti_sec 
        RR_pp1(r,:) = [];
    end
end 
fine1 = RR_pp1(end,2)-cinqueminuti_sec;
conta1 = 0;
for r = length(RR_pp1):-1:1
    if RR_pp1(r,2)>fine1
        conta1 = conta1 + 1;
        RR_pp1(r,:) = [];
    end 
end

figure()
plot(RR(:,2),RR(:,1)); hold on;plot(RR_pp1(:,2),RR_pp1(:,1),'r'),xlabel('[s]'),ylabel('RR[s]'),legend('Pre Processing','Post Processing'),...
    title('RR interval - without first and last 5 minutes')
figure(), subplot(211),plot(RR_pp1(:,2),RR_pp1(:,1)), subplot(212),boxplot(RR_pp1(:,1)),sgtitle('RR interval - without first and last 5 minutes')
        
figure()
subplot(121), plot(RR_tempo,RR_sec); hold on;plot(RR_pp(:,2),RR_pp(:,1),'r'),xlabel('[s]'),ylabel('RR[s]'),legend('Pre Processing','Post Processing'),...
    title('RR interval')
subplot(122), plot(RR(:,2),RR(:,1)); hold on;plot(RR_pp1(:,2),RR_pp1(:,1),'r'),xlabel('[s]'),ylabel('RR[s]'),legend('Pre Processing','Post Processing'),...
    title('RR interval - without first and last 5 minutes')