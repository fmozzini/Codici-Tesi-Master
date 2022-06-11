%% SCG Peaks - Z axis filtered acceleration
%% RICERCA MODIFICATA: uso una finestra di 10 s presa durante la notte
clear all
close all
clc

%%
folderSCG = 'C:\Users\feder\Desktop\Tesi\Data\Filtered Signals';

list = dir(folderSCG);
list(1) = [];
list(1) = [];
N = length(list);
list(N-1) = [];
list(N-1) = [];
% N = length(list)
N = length(list)-1;

[~,txtdata] = xlsread('C:\Users\feder\Desktop\Tesi\Info Pazienti.xlsx','H:I');
txtdata(1,:) = [];
Inizio_Holter = txtdata(:,1);
Periodo_Sonno = txtdata(:,2);
addpath 'C:\Users\feder\Desktop\Tesi'\Data\'Filtered Signals'\
% addpath 'C:\Users\feder\Desktop\Tesi'\Data\'PostProc SCG'\
addpath 'C:\Users\feder\Desktop\Tesi'\Codes\
%% 
for m = 22:N

    FOLDERSCG = fullfile(list(m).folder, list(m).name)
    file = dir(FOLDERSCG);
    name = file.name;
    load(name)


    dim = size(Acc_filt.z_filt,1); %durata in campioni
    fs_SCG = 64;
    wind_30sec = fs_SCG*30;
    n_window = round(dim/wind_30sec);
    % template = [];
    HR_picchi = zeros(dim,1);
    Pos_picchi = zeros(dim,1);
    Amp_picchi = zeros(dim,1);
    n_peaks = 0;
    peaksFORwindow_SCG30 = zeros(n_window,1); 

    % Template
    % Cerco la notte 
    Hour_add = 2*3600;
    Inizio = datevec(Inizio_Holter{m}); 
    Sonno = datevec(Periodo_Sonno{m}); 
    durata_sec = etime(Sonno,Inizio)
    durata_h = durata_sec/3600
    Night = durata_h + 2;
    Night_sec = round(Night*3600*64);
%     templateNight = Acc_filt.z_filt(Night_sec:Night_sec+30*fs_SCG);
%     plot(Night_sec:Night_sec+30*fs_SCG,templateNight),title('Template di 30 sec preso durante la notte')
    for n = 1:2:200
        templateNight = Acc_filt.z_filt(Night_sec+10*n*fs_SCG:Night_sec+10*(n+1)*fs_SCG);
        plot(Night_sec+10*n*fs_SCG:Night_sec+10*(n+1)*fs_SCG,templateNight),title('Template di 30 sec preso durante la notte')
        pause
        answer = questdlg('Do you want to use this template?', ...
	        'SCG template matching template ', ...
	        'Yes','No','Yes');
        switch answer
            case 'No'
                continue
            case 'Yes'
                TemplateNight = templateNight;
                break 
        end 
    end 
%     [~,~,~,Template] = SCG_template_matching_corr(TemplateNight,fs_SCG,'x',[]);
%     plot(Template)
%     pause
    for i = 1:n_window-1 
        [~, pos_picchi, amp_picchi,~] = SCG_template_matching_corr(Acc_filt.z_filt((i-1)*wind_30sec+1:wind_30sec*i),fs_SCG,'x',TemplateNight);
        n_peaks_new = size(pos_picchi,1);
        if i == 1
            Pos_picchi(n_peaks+1:n_peaks+n_peaks_new,:) = pos_picchi;
        else
            Pos_picchi(n_peaks+1:n_peaks+n_peaks_new,:) = pos_picchi+(i-1)*wind_30sec ; 
        end 
        Amp_picchi(n_peaks+1:n_peaks+n_peaks_new,:) = amp_picchi;
        n_peaks = n_peaks + n_peaks_new;
        peaksFORwindow_SCG30(i) = n_peaks;  
    end 
    num = n_peaks;
    AMP_picchi_SCG = zeros(num,1);
    POS_picchi_SCG = zeros(num,1);
    AMP_picchi_SCG = Amp_picchi(1:num);
    POS_picchi_SCG = Pos_picchi(1:num); 
    
    HR_SCG = 60./(diff(POS_picchi_SCG))*fs_SCG;

    % Save
    name_SCG = erase(name,"FILT-")
    Acc_z = Acc_filt.z_filt;
    save(['C:\Users\feder\Desktop\Tesi\Data\Picchi SCG - Acc z - 10 sec\' 'SCG(Az)_picchi_10SEC-' name_SCG],'Acc_z','AMP_picchi_SCG','POS_picchi_SCG','peaksFORwindow_SCG30','HR_SCG')
    save(['C:\Users\feder\Desktop\Tesi\Data\Picchi SCG - Acc z - 10 sec\Template\' 'template_10SEC-' name_SCG],'templateNight')
end 