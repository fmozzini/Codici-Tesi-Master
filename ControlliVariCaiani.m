clear all
close all
clc

%%
folderWID = 'C:\Users\feder\Desktop\Tesi\Data\Windows_10SEC';
folderPT = 'C:\Users\feder\Desktop\Tesi\Data\PostProc PT_10SEC'; 
folderECG = 'C:\Users\feder\Desktop\Tesi\Data\Filtered ECG';
folderPAR = 'C:\Users\feder\Desktop\Tesi\Data\Parameters SCG_10SEC';
folderSCG = 'C:\Users\feder\Desktop\Tesi\Data\Picchi SCG - Acc z - 10 sec';
folderFP = 'C:\Users\feder\Desktop\Tesi\Data\Fiducial Points SCG_10SEC';

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
listFP = dir(folderFP);
listFP(1) = [];
listFP(1) = [];

addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Filtered ECG'\
addpath 'C:\Users\feder\Desktop\Tesi'\Data\'PostProc PT_10SEC'\
addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Windows_10SEC'\
addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Parameters SCG_10SEC'\
addpath 'C:\Users\feder\Desktop\Tesi'\Codes\
addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Picchi SCG - Acc z - 10 sec'\
addpath 'C:\Users\feder\Desktop\Tesi\Data\Fiducial Points SCG_10SEC\';
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

    FOLDERFP = fullfile(listFP(m).folder, listFP(m).name)
    file = dir(FOLDERFP);
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
     R0 = R(pos_R0,:);MC0 = MC(pos_R0,:); IVC0 = IVC(pos_R0,:); AO0 = AO(pos_R0,:); minAO_RE0 = minAO_RE(pos_R0,:); RE0 = RE(pos_R0,:); AC0 = AC(pos_R0,:); minbeforeAC0 = minbeforeAC(pos_R0,:); 
     T0 = T(pos_R0,:); fine_T0 = fine_T(pos_R0,:); Q0 = Q(pos_R0,:);
     pos_R4 = find(R(:,3)==4);
     R4 = R(pos_R4,:);MC4 = MC(pos_R4,:); IVC4 = IVC(pos_R4,:); AO4 = AO(pos_R4,:); minAO_RE4 = minAO_RE(pos_R4,:); RE4 = RE(pos_R4,:); AC4 = AC(pos_R4,:); minbeforeAC4 = minbeforeAC(pos_R4,:); 
     T4 = T(pos_R4,:); fine_T4 = fine_T(pos_R4,:); Q4 = Q(pos_R4,:);
     pos_R5 = find(R(:,3)==5);
     R5 = R(pos_R5,:);MC5 = MC(pos_R5,:); IVC5 = IVC(pos_R5,:); AO5 = AO(pos_R5,:); minAO_RE5 = minAO_RE(pos_R5,:); RE5 = RE(pos_R5,:); AC5 = AC(pos_R5,:); minbeforeAC5 = minbeforeAC(pos_R5,:); 
     T5 = T(pos_R5,:); fine_T5 = fine_T(pos_R5,:); Q5 = Q(pos_R5,:);
% PER ESSERE SICURO DELLE RIGHE CHE PRENDO, DEVO SCRIVERMI DA QUALCHE PARTE
% A CHE RIGA CORRISPONDEEEEEEEEEEEE
     N = 100;
     randomIndexes0 = randperm(length(R0), N);
     randomIndexes4 = randperm(length(R4), N);
     randomIndexes5 = randperm(length(R5), N);
     R0_random = R0(randomIndexes0,:); MC0_random = MC0(randomIndexes0,:); IVC0_random = IVC0(randomIndexes0,:); AO0_random = AO0(randomIndexes0,:); minAO_RE0_random = minAO_RE0(randomIndexes0,:);
     RE0_random = RE0(randomIndexes0,:); AC0_random = AC0(randomIndexes0,:); minbeforeAC0_random = minbeforeAC0(randomIndexes0,:); T0_random = T0(randomIndexes0,:); fine_T0_random = fine_T0(randomIndexes0,:);
     Q0_random = Q0(randomIndexes0,:);
     R4_random = R4(randomIndexes4,:); MC4_random = MC4(randomIndexes4,:); IVC4_random = IVC4(randomIndexes4,:); AO4_random = AO4(randomIndexes4,:); minAO_RE4_random = minAO_RE4(randomIndexes4,:);
     RE4_random = RE4(randomIndexes4,:); AC4_random = AC4(randomIndexes4,:); minbeforeAC4_random = minbeforeAC4(randomIndexes4,:); T4_random = T4(randomIndexes4,:); fine_T4_random = fine_T4(randomIndexes4,:);
     Q4_random = Q4(randomIndexes4,:);
    
     R5_random = R5(randomIndexes5,:); MC5_random = MC5(randomIndexes5,:); IVC5_random = IVC5(randomIndexes5,:); AO5_random = AO5(randomIndexes5,:); minAO_RE5_random = minAO_RE5(randomIndexes5,:);
     RE5_random = RE5(randomIndexes5,:); AC5_random = AC5(randomIndexes5,:); minbeforeAC5_random = minbeforeAC5(randomIndexes5,:); T5_random = T5(randomIndexes5,:); fine_T5_random = fine_T5(randomIndexes5,:);
     Q5_random = Q5(randomIndexes5,:);
%      R_random = [ R0_random; R4_random; R5_random]; MC_random = [ MC0_random; MC4_random; MC5_random]; IVC_random = [ IVC0_random; IVC4_random; IVC5_random]; AO_random = [ AO0_random; AO4_random; AO5_random];
%      minAO_RE_random = [ minAO_RE0_random; minAO_RE4_random; minAO_RE5_random]; RE_random = [ RE0_random; RE4_random; RE5_random]; AC_random = [ AC0_random; AC4_random; AC5_random];
%      minbeforeAC_random = [ minbeforeAC0_random; minbeforeAC4_random; minbeforeAC5_random]; T_random = [ T0_random; T4_random; T5_random]; fine_T_random = [ fine_T0_random; fine_T4_random; fine_T5_random];
%      Q_random = [ Q0_random; Q4_random; Q5_random];
%      R_random1 = sortrows(R_random); MC_random1 = sortrows(MC_random); IVC_random1 = sortrows(IVC_random); AO_random1 = sortrows(AO_random); minAO_RE_random1 = sortrows(minAO_RE_random); RE_random1 = sortrows(RE_random);
%      AC_random1 = sortrows(AC_random); minbeforeAC_random1 = sortrows(minbeforeAC_random); T_random1 = sortrows(T_random); fine_T_random1 = sortrows(fine_T_random); Q_random1 = sortrows(Q_random);   
    R_random1 = sortrows(R5_random); MC_random1 = sortrows(MC5_random); IVC_random1 = sortrows(IVC5_random); AO_random1 = sortrows(AO5_random); minAO_RE_random1 = sortrows(minAO_RE5_random); RE_random1 = sortrows(RE5_random);
     AC_random1 = sortrows(AC0_random); minbeforeAC_random1 = sortrows(minbeforeAC5_random); T_random1 = sortrows(T5_random); fine_T_random1 = sortrows(fine_T5_random); Q_random1 = sortrows(Q5_random);   
     for i = 1:length(R_random1)
%             count = count+1;
            R_succ = R(find(R==R_random1(i))+1,:);
            qrs1 = (R_random1(i,1)/1024)*64;
            qrs2 = (R_succ(1)/1024)*64;
            finestrabattito_ECGF = ECG_filt(R_random1(i,1)-window_ECG:R_succ(1)-window_ECG)';
            finestrabattito_SCG = Acc_z(qrs1-window_SCG:qrs2-window_SCG)';
            figure()
            set(gcf, 'WindowState', 'maximized');
            a = subplot(211); plot((R_random1(i,1)-window_ECG:R_succ(1)-window_ECG)./1024,finestrabattito_ECGF),xlabel('[s]'),title('Ecg'); hold on;
%             plot(R_random1(i,1)./1024,R_random1(i,2),'*r'); hold on; plot(T_random1(i,1)/1024,T_random1(i,2),'*g'); hold on;
%             xline(R_random1(i,1)/1024); hold on; xline(T_random1(i,1)/1024,'--g'); hold on; plot(Q_random1(i,1)/1024,Q_random1(i,2),'*b'); hold on;
%             plot(fine_T_random1(i,1)/1024,fine_T_random1(i,2),'*b'); hold on; xline(fine_T_random1(i,1)/1024,'--b')
%             text(R_random1(i,1)./1024,R_random1(i,2),' R')
%             text(T_random1(i,1)/1024,T_random1(i,2),' T')
%             text(fine_T_random1(1)/1024,fine_T_random1(2),' fine onda T')
%             text(Q_random1(i,1)/1024,Q_random1(i,2),' Q')
            b = subplot(212); plot((qrs1-window_SCG:qrs2-window_SCG)./64,finestrabattito_SCG),xlabel('[s]'),title('SCG'); hold on;
%             xline(qrs1/64); hold on; xline(T_random1(i,1)/1024,'--g'); hold on; 
%             xline(fine_T_random1(i,1)/1024,'--b'); hold on;
%             xline((AO_random1(i,1)+window_SCG)./64); hold on 
%             plot(IVC_random1(i,1)/64,IVC_random1(i,2),'*r'); hold on
%             plot(AO_random1(i,1)/64,AO_random1(i,2),'*r'); hold on;
%             plot(RE_random1(i,1)/64,RE_random1(i,2),'*r'); hold on;
%             plot(MC_random1(i,1)/64,MC_random1(i,2),'*r'); hold on;
% %             line([AO_random1(i,1)/64 IVC_random1(i,1)/64],[AO_random1(i,2) IVC_random1(i,2)],'Color','red','LineStyle','--'); hold on;
%             plot(minAO_RE_random1(i,1)/64,minAO_RE_random1(i,2),'*r'); hold on;
%             text(IVC_random1(i,1)/64,IVC_random1(i,2),' IVC')
%             text(AO_random1(i,1)/64,AO_random1(i,2),' AO')
%             text(RE_random1(i,1)/64,RE_random1(i,2),' RE')
%             text(MC_random1(i,1)/64,MC_random1(i,2),' MC')
%             text(minAO_RE_random1(i,1)/64,minAO_RE_random1(i,2),' min AO-RE');
%             if R_random1(i,3) == 0
%                 plot(AC_random1(i,1)/64,AC_random1(i,2),'*r'); hold on;
%                 plot(minbeforeAC_random1(i,1)/64,minbeforeAC_random1(i,2),'*r');
%                 text(AC_random1(i,1)/64,AC_random1(i,2),' AC')
%                 text(minbeforeAC_random1(i,1)/64,minbeforeAC_random1(i,2), 'min before AC')
%             end 
            sgtitle(find(R==R_random1(i)))
            pause
            sgtitle(R_random1(i,3))
            pause
            close all
    end 
  




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

