%% Cut signals - Elimino i primi 5 minuti dei segnali, sono troppo corrotti 
 clear all
 close all
 clc

folderECG = 'C:\Users\feder\Desktop\Tesi\Data\Filtered ECG';
folderSCG = 'C:\Users\feder\Desktop\Tesi\Data\Filtered Signals';
listECG = dir(folderECG);
listECG(1) = [];
listECG(1) = [];
listECG(end) = [];
list = dir(folderSCG);
list(1) = [];
list(1) = [];
N = length(list);
list(N-1) = [];
list(N-1) = [];
N = length(list)

addpath 'C:\Users\feder\Desktop\Tesi'\Data\''\'Filtered ECG'\
addpath 'C:\Users\feder\Desktop\Tesi\Data\Filtered Signals'\
addpath 'C:\Users\feder\Desktop\Tesi'\Codes\

%% 
for m = 1:1
    FOLDERSCG = fullfile(list(m).folder, list(m).name)
    file = dir(FOLDERSCG);
    name = file.name;
    load(name)

    FOLDERECG = fullfile(listECG(m).folder, listECG(m).name)
    file = dir(FOLDERECG);
    name = file.name;
    load(name)

    fs_ECG = 1024;
    fs_SCG = 64;
    
    cinqueminuti = 5*60 %300 sec 
    cinqueminuti_ECG = cinqueminuti*fs_ECG; %307200
    cinqueminuti_SCG = cinqueminuti*fs_SCG; %19200


    ECG_filt = ECG_filt(cinqueminuti_ECG+1:length(ECG_filt));
    Acc_x = Acc_filt.x_filt(cinqueminuti_SCG+1:length(Acc_filt.x_filt));
    Acc_y = Acc_filt.y_filt(cinqueminuti_SCG+1:length(Acc_filt.y_filt));
    Acc_z = Acc_filt.z_filt(cinqueminuti_SCG+1:length(Acc_filt.z_filt));
    Acc_tempo = Acc_filt.Var4(cinqueminuti_SCG+1:length(Acc_filt.Var4));
    Rot_x = Rot_filt.x_rotfilt(cinqueminuti_SCG+1:length(Rot_filt.x_rotfilt));
    Rot_y = Rot_filt.y_rotfilt(cinqueminuti_SCG+1:length(Rot_filt.y_rotfilt));
    Rot_z = Rot_filt.z_rotfilt(cinqueminuti_SCG+1:length(Rot_filt.z_rotfilt));
    Rot_tempo = Rot_filt.Var4(cinqueminuti_SCG+1:length(Rot_filt.Var4));
    
    Acc_filt = table(Acc_x,Acc_y,Acc_z,Acc_tempo);
    Acc_filt.Properties.VariableNames = {'x_filt' 'y_filt' 'z_filt' 'Var4'};
    Rot_filt = table(Rot_x,Rot_y,Rot_z,Rot_tempo);
    Rot_filt.Properties.VariableNames = {'x_rotfilt' 'y_rotfilt' 'z_rotfilt' 'Var4'};

    name = erase(name,"ECG_FILT-")
    save(['C:\Users\feder\Desktop\Tesi\Data\PostProc SCG\' 'PP_SCG -' name],'Acc_filt','Rot_filt')
    save(['C:\Users\feder\Desktop\Tesi\Data\PostProc ECG\' 'PP_ECG -' name],'ECG_filt')
end 