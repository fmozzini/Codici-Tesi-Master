clear all
close all
clc
%%
folderFP = 'C:\Users\feder\Desktop\Tesi\Data\Fiducial Points SCG';
folderECG = 'C:\Users\feder\Desktop\Tesi\Data\Filtered ECG';
folderPT = 'C:\Users\feder\Desktop\Tesi\Data\PostProc PT 1'; 

listPT = dir(folderPT);
listPT(1) = [];
listPT(1) = [];
listECG = dir(folderECG);
listECG(1) = [];
listECG(1) = [];
listECG(end) = [];
listFP = dir(folderFP);
listFP(1) = [];
listFP(1) = [];
N = length(listPT);

addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Filtered ECG'\
addpath 'C:\Users\feder\Desktop\Tesi'\Data\'PostProc PT 1'\
addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Fiducial Points SCG'\
addpath 'C:\Users\feder\Desktop\Tesi'\Codes\
%%
for m = 1:1
    FOLDERFP = fullfile(listFP(m).folder, listFP(m).name)
    file = dir(FOLDERFP);
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

    for i = 1: % devo scorrere tutti i picchi! POSSO METTERE UNA CONDIZIONE SU NUMERO DI PICCHI/GIORNO E NOTTE 

    figure()
    a = subplot(211)
    plot((qrs_I(i)-window_ECG:qrs_I(i+1)-window_ECG)./1024,finestrabattito_ECG),xlabel('[s]'); hold on; plot(qrs_I(i)./1024,qrs_AMP(i),'*r'); hold on; plot(T(i,1)/1024,T(i,2),'*g'); hold on;
    xline(qrs_I(i)/1024); hold on; xline(T(i,1)/1024,'--g'); hold on;
    xline(T40d/1024,'--m'); hold on; xline(T80d/1024,'--m'); hold on;
    plot(Xmd/1024,ym,'*m'); hold on; plot(Xrd/1024,yr,'*m'); hold on;
    plot(Xid/1024,yi,'*b'); hold on;
    plot(Q(i,1)/1024,Q(i,2),'*b'); hold on;
%             xline(QS2maxline_ECGd/1024,'--r'); hold on;
    xline(RR_sysd/1024,'--y'); hold on;
    plot(fine_T(i,1)/1024,fine_T(i,2),'*b');
    hold on; xline(fine_T(i,1)/1024,'--b')
    text(qrs_I(i)./1024,qrs_AMP(i),' R')
    text(T(i,1)/1024,T(i,2),' T')
    text(fine_T(1)/1024,fine_T(2),' fine onda T')
    text(Q(i,1)/1024,Q(i,2),' Q')
    b = subplot(212)
    plot((qrs1-window_SCG:qrs2-window_SCG)./64,finestrabattito_SCG),xlabel('[s]'); hold on; 
    for r = 1:length(row)
        plot((POS_picchi_SCG(row(r)))./64,AMP_picchi_SCG(row(r)),'mo')
    end 
    xline(qrs1/64); hold on; xline(T(i,1)/1024,'--g'); hold on; 
    xline(fine_T(i,1)/1024,'--b'); hold on;
    xline((AO(i,1)+window_SCG)./64); hold on 
    xline(QS2maxlined/64,'--m'); hold on;
    xline((qrs1+0.005*64)/64,'r'); hold on;
    plot(IVC(i,1)/64,IVC(i,2),'*r'); hold on
    plot(AO(i,1)/64,AO(i,2),'*r'); hold on;
    plot(RE(i,1)/64,RE(i,2),'*r'); hold on;
    plot(AC(i,1)/64,AC(i,2),'*r'); hold on;
    plot(MC(i,1)/64,MC(i,2),'*r'); hold on;
    line([AO(i,1)/64 IVC(i,1)/64],[AO(i,2) IVC(i,2)],'Color','red','LineStyle','--'); hold on;
    plot(minAO_RE(i,1)/64,minAO_RE(i,2),'*r'); hold on;
    plot(minbeforeAC(i,1)/64,minbeforeAC(i,2),'*r');
    text(IVC(i,1)/64,IVC(i,2),' IVC')
    text(AO(i,1)/64,AO(i,2),' AO')
    text(RE(i,1)/64,RE(i,2),' RE')
    text(AC(i,1)/64,AC(i,2),' AC')
    text(MC(i,1)/64,MC(i,2),' MC')
    text(minAO_RE(i,1)/64,minAO_RE(i,2),' min AO-RE');
    text(minbeforeAC(i,1)/64,minbeforeAC(i,2), 'min before AC')
    sgtitle(i)
    R_MC1(i,1)
     pause
     close all
end 