clear all
close all
clc

%%
folderWID = 'C:\Users\feder\Desktop\Tesi\Data\Windows';
folderPT = 'C:\Users\feder\Desktop\Tesi\Data\PostProc PT'; 
folderECG = 'C:\Users\feder\Desktop\Tesi\Data\Filtered ECG';
folderPAR = 'C:\Users\feder\Desktop\Tesi\Data\Parameters SCG';
folderSCG = 'C:\Users\feder\Desktop\Tesi\Data\Picchi SCG - Acc z';

listWIND = dir(folderWID);
listWIND(1) = [];
listWIND(1) = [];
listPT = dir(folderPT);
listPT(1) = [];
listPT(1) = [];
listECG = dir(folderECG);
listECG(1) = [];
listECG(1) = [];
listECG(end) = [];
listPAR = dir(folderPAR);
listPAR(1) = [];
listPAR(1) = [];
listSCG = dir(folderSCG);
listSCG(1) = [];
listSCG(1) = [];


addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Filtered ECG'\
addpath 'C:\Users\feder\Desktop\Tesi'\Data\'PostProc PT'\
addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Windows'\
addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Parameters SCG'\
addpath 'C:\Users\feder\Desktop\Tesi'\Codes\
addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Picchi SCG - Acc z'\
%%
for m = 1:1
    FOLDERWIND = fullfile(listWIND(m).folder, listWIND(m).name)
    file = dir(FOLDERWIND);
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

    FOLDERPAR = fullfile(listPAR(m).folder, listPAR(m).name)
    file = dir(FOLDERPAR);
    name = file.name;
    load(name)

    FOLDERSCG = fullfile(listSCG(m).folder, listSCG(m).name)
    file = dir(FOLDERSCG);
    name = file.name;
    load(name)

%% Conto quanti picchi per tipologia
    n_picchi0 = 0;
    n_picchi1 = 0;
    n_picchi2 = 0;
    n_picchi3 = 0;
    n_picchipiudi3 = 0;
    for i = 1:length(n_picchi)
        if n_picchi(i) == 0
            n_picchi0 = n_picchi0 + 1;
        elseif  n_picchi(i) == 1
            n_picchi1 = n_picchi1 + 1;
        elseif n_picchi(i) == 2
            n_picchi2 = n_picchi2 + 1;
        elseif n_picchi(i) == 3
            n_picchi3 = n_picchi3 + 1;     
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
% GN = 1 notte; GN = 0 giorno
    fs_ECG = 1024;
    fs_SCG = 64;
    window_SCG = 0.2*fs_SCG;
    window_ECG = 0.2*fs_ECG;
    count = 0;
    for i = 1:length(R)-1
        if n_picchi(i) == 0 && GN(i) == 1
            count = count+1;
            qrs1 = (R(i,1)/1024)*64;
            qrs2 = (R(i+1,1)/1024)*64;
            finestra_dopo = qrs2-window_SCG; %cambiato
            finestra_prima = qrs1-window_SCG;
            [row]=find(POS_picchi_SCG<finestra_dopo & POS_picchi_SCG>finestra_prima);
            finestrabattito_ECGF = ECG_filt(R(i,1)-window_ECG:R(i+1,1)-window_ECG)';
            finestrabattito_SCG = Acc_z(qrs1-window_SCG:qrs2-window_SCG)';
            figure()
            set(gcf, 'WindowState', 'maximized');
            a = subplot(211); plot((R(i,1)-window_ECG:R(i+1,1)-window_ECG)./1024,finestrabattito_ECGF),xlabel('[s]'),title('Ecg'); hold on;
            plot(R(i,1)/1024,R(i,2),'*r')
            b = subplot(212); plot((qrs1-window_SCG:qrs2-window_SCG)./64,finestrabattito_SCG),xlabel('[s]'),title('SCG'); hold on;
            for r = 1:length(row)
                plot((POS_picchi_SCG(row(r)))./64,AMP_picchi_SCG(row(r)),'mo')
            end
            sgtitle(i)
            pause
            close all
        end 
    end 

    %% RIESCO A RECUPERARE I PICCHI PT IN FINESTRE CON PIU' DI 2 PICCHI (>= 3)
    % GN = 1 notte; GN = 0 giorno
    fs_ECG = 1024;
    fs_SCG = 64;
    window_SCG = 0.2*fs_SCG;
    window_ECG = 0.2*fs_ECG;
    count = 0;
    qrs_new = [];
    error = 0;
    for i = 1:length(R)-1
        if n_picchi(i) >= 3 
            count = count+1;
            qrs1 = (R(i,1)/1024)*64;
            qrs2 = (R(i+1,1)/1024)*64;
            finestra_dopo = qrs2-window_SCG; %cambiato
            finestra_prima = qrs1-window_SCG;
            [row]=find(POS_picchi_SCG<finestra_dopo & POS_picchi_SCG>finestra_prima);
            finestrabattito_ECGF = ECG_filt(R(i,1)-window_ECG:R(i+1,1)-window_ECG)';
            finestrabattito_SCG = Acc_z(qrs1-window_SCG:qrs2-window_SCG)';
            try 
            [~,qrs_i_raw,delay]=pan_tompkin(finestrabattito_ECGF,fs_ECG,0);
            n_peaks_new = size(qrs_i_raw,2);
            inizio = R(i,1)-qrs_i_raw(1);
            qrs_i_new = inizio+qrs_i_raw;
            qrs_amp_new = ECG_filt(qrs_i_new);
            qrs_new1 = [qrs_i_new' qrs_amp_new];
            qrs_new = [qrs_new; qrs_new1];
%             figure()
%             set(gcf, 'WindowState', 'maximized');
%             a = subplot(211); plot((R(i,1)-window_ECG:R(i+1,1)-window_ECG)./1024,finestrabattito_ECGF),xlabel('[s]'),title('Ecg'); hold on;
%             plot(R(i,1)/1024,R(i,2),'*r'); hold on;
%             for p = 1: length(qrs_i_new)
%                 plot(qrs_new1(p,1)/1024,qrs_new1(p,2),'*g')
%             end 
%             b = subplot(212); plot((qrs1-window_SCG:qrs2-window_SCG)./64,finestrabattito_SCG),xlabel('[s]'),title('SCG'); hold on;
%             for r = 1:length(row)
%                 plot((POS_picchi_SCG(row(r)))./64,AMP_picchi_SCG(row(r)),'mo')
%             end
%             sgtitle(i)
%             pause 
%             close all
            catch ME 
                error = error+1;
            end
        end 
    end 
end 

 %% Select random 100 indixes for each tag and plot 
     fs_ECG = 1024;
     fs_SCG = 64;
     window_SCG = 0.2*fs_SCG;
     window_ECG = 0.2*fs_ECG;
     pos_R0 = find(R(:,3)==0);
     R0 = R(pos_R0,:);
     pos_R4 = find(R(:,3)==4);
     R4 = R(pos_R4,:);
     pos_R5 = find(R(:,3)==5);
     R5 = R(pos_R5,:);
     N = 100;
     randomIndexes0 = randperm(length(R0), N);
     randomIndexes4 = randperm(length(R4), N);
     randomIndexes5 = randperm(length(R5), N);
     R0_random = R0(randomIndexes0,:);
     R4_random = R4(randomIndexes4,:);
     R5_random = R5(randomIndexes5,:);
     R_random = [ R0_random; R4_random; R5_random];
     R_random1 = sortrows(R_random);
     for i = 1:length(R_random1)
            count = count+1;
            qrs1 = (R_random1(i,1)/1024)*64;
            qrs2 = (R_random1(i+1,1)/1024)*64;
            finestrabattito_ECGF = ECG_filt(R_random1(i,1)-window_ECG:R_random1(i+1,1)-window_ECG)';
            finestrabattito_SCG = Acc_z(qrs1-window_SCG:qrs2-window_SCG)';
            figure()
            set(gcf, 'WindowState', 'maximized');
            a = subplot(211); plot((R_random1(i,1)-window_ECG:R_random1(i+1,1)-window_ECG)./1024,finestrabattito_ECGF),xlabel('[s]'),title('Ecg')
            b = subplot(212); plot((qrs1-window_SCG:qrs2-window_SCG)./64,finestrabattito_SCG),xlabel('[s]'),title('SCG')
            sgtitle(R_random1(i,3))
            pause
            close all
    end 
  


% %% Conto il numero di righe consecutive per tag0 e tag 5 (quanti battiti consecutivi ho durante il giorno e durante la notte)
% % Tag 0, Notte
% Maxduratanotte_0 = 0;
% duratanotte_0 = 0;
% for i = 2:length(R)-1
%     if R(i,3) == 0 && GN(i) == 0 %notte
%         if R(i-1,3) == 0 && GN(i-1) == 0
%             duratanotte_0 = duratanotte_0+1;
%             if duratanotte_0 > Maxduratanotte_0
%                 Maxduratanotte_0 = duratanotte_0;
%             end 
%         else
%             duratanotte_0 = 0;
%         end 
%     end
% end 
% 
% % Tag 0, Giorno
% Maxduratagiorno_0 = 0;
% duratagiorno_0 = 0;
% for i = 2:length(R)-1
%     if R(i,3) == 0 && GN(i) == 1 %giorno
%         if R(i-1,3) == 0 && GN(i-1) == 1
%             duratagiorno_0 = duratagiorno_0+1;
%             if duratagiorno_0 > Maxduratagiorno_0
%                 Maxduratagiorno_0 = duratagiorno_0;
%             end 
%         else
%             duratagiorno_0 = 0;
%         end 
%     end
% end 
% %%
% % Tag 0 e 5, Giorno
% Maxduratagiorno_05 = 0;
% duratagiorno_05 = 0;
% for i = 2:length(R)-1
%     if R(i,3) == 0 || R(i,3) == 5
%         if GN(i) == 1 %giorno
%             if R(i-1,3) == 0 || R(i-1,3) == 5
%                 if GN(i-1) == 1
%                     duratagiorno_05 = duratagiorno_05+1;
%                     if duratagiorno_05 > Maxduratagiorno_05
%                         Maxduratagiorno_05 = duratagiorno_05;
%                     end 
%                 else
%                     duratagiorno_05 = 0;
%                 end 
%             end 
%         end 
%     end
% end 
% 
% % Tag 0 e 5, Notte
% Maxduratanotte_05 = 0;
% duratanotte_05 = 0;
% for i = 2:length(R)-1
%     if R(i,3) == 0 || R(i,3) == 5
%         if GN(i) == 0 %notte
%             if R(i-1,3) == 0 || R(i-1,3) == 5
%                 if GN(i-1) == 0
%                     duratanotte_05 = duratanotte_05+1;
%                     if duratanotte_05 > Maxduratanotte_05
%                         Maxduratanotte_05 = duratanotte_05;
%                     end 
%                 else
%                     duratanotte_05 = 0;
%                 end 
%             end 
%         end 
%     end
% end 
% %% FUNZIAAAA - Ema Gero, aiutante Crocco 
% M = zeros(length(R), 2);
% c = 0;
% ln = 1;
% val = R(1,3);
% for i = 2:length(R)
% 
%     c = c + 1;
%     
%     cur_val = R(i,3);
%     if cur_val ~= val
%         M(ln, 1) = val;
%         M(ln, 2) = c;
%         c = 0;
%         val = cur_val;
%         ln = ln+1;
%     end
% 
% end
% 
% M(ln, 1) = val;
% M(ln, 2) = c;
% M(ln:end,:) = [];

%%  
    load 'template-2021-01-16 15.01.53.mat'
    fs_ECG = 1024;
    fs_SCG = 64;
    window_SCG = 0.2*fs_SCG;
    window_ECG = 0.2*fs_ECG;
    count = 0;
    for i = 1:length(R)
        if n_picchi(i) >=3 1 && GN(i) == 1
            count = count+1;
            qrs1 = (R(i,1)/1024)*64;
            cinquesecECG = 5*1024;
            cinquesecSCG = 5*64;
            finestrabattito_ECGF = ECG_filt(R(i,1)-cinquesecECG:R(i,1)+cinquesecECG)';
            finestrabattito_SCG = Acc_z(qrs1-cinquesecSCG:qrs1+cinquesecSCG)';
            figure()
            set(gcf, 'WindowState', 'maximized');
            a = subplot(211); plot((R(i,1)-cinquesecECG:R(i,1)+cinquesecECG)./1024,finestrabattito_ECGF),xlabel('[s]'),title('Ecg'); hold on;
%             plot(R(i,1)/1024,R(i,2),'*r')
            b = subplot(212); plot((qrs1-cinquesecSCG:qrs1+cinquesecSCG)./64,finestrabattito_SCG),xlabel('[s]'),title('SCG'); hold on;
%             for r = 1:length(row)
%                 plot((POS_picchi_SCG(row(r)))./64,AMP_picchi_SCG(row(r)),'mo')
%             end
            hold on; plot((qrs1-cinquesecSCG:qrs1+cinquesecSCG)./64,template,'r')
            sgtitle(i)
            pause
            close all
        end 
    end 

%% Provo a plottare tutti i template che uso durante SCG 
clear all
close all
clc

%%
folderT = 'C:\Users\feder\Desktop\Tesi\Data\Picchi SCG - Acc z\Template';
listT = dir(folderT);
listT(1) = [];
listT(1) = [];
N = length(listT);
TEMPLATE = zeros(N,641);
for m = 1:N
    FOLDERT = fullfile(listT(m).folder, listT(m).name)
    file = dir(FOLDERT);
    name = file.name;
    load(name)
    TEMPLATE(m,:) = template';
end 
figure()
subplot(331),plot(TEMPLATE(1,:)),title('1'); subplot(332),plot(TEMPLATE(2,:)),title('2'); subplot(333),plot(TEMPLATE(3,:)),title('3'); 
subplot(334),plot(TEMPLATE(4,:)),title('4'); subplot(335),plot(TEMPLATE(5,:)),title('5'); subplot(336),plot(TEMPLATE(6,:)),title('6');
subplot(337),plot(TEMPLATE(7,:)),title('7'); subplot(338),plot(TEMPLATE(8,:)),title('8'); subplot(339),plot(TEMPLATE(9,:)),title('9'); 

figure()
subplot(331),plot(TEMPLATE(10,:)),title('10'); subplot(332),plot(TEMPLATE(11,:)),title('11'); subplot(333),plot(TEMPLATE(12,:)),title('12'); 
subplot(334),plot(TEMPLATE(13,:)),title('13'); subplot(335),plot(TEMPLATE(14,:)),title('14'); subplot(336),plot(TEMPLATE(15,:)),title('15');
subplot(337),plot(TEMPLATE(16,:)),title('16'); subplot(338),plot(TEMPLATE(17,:)),title('17'); subplot(339),plot(TEMPLATE(18,:)),title('18'); 

figure()
subplot(331),plot(TEMPLATE(19,:)),title('19'); subplot(332),plot(TEMPLATE(20,:)),title('20'); subplot(333),plot(TEMPLATE(21,:)),title('21'); 
subplot(334),plot(TEMPLATE(22,:)),title('22'); subplot(335),plot(TEMPLATE(23,:)),title('23'); subplot(336),plot(TEMPLATE(24,:)),title('24');
subplot(337),plot(TEMPLATE(25,:)),title('25'); subplot(338),plot(TEMPLATE(26,:)),title('26');

