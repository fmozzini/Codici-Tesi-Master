% FilterDATA    
clear all
close all
clc

%% 
folder = 'C:\Users\feder\Desktop\Tesi\Data\24h Signals';
list = dir(folder);
list(1) = [];
list(1) = [];
N = length(list);
list(N-1) = [];
list(N-1) = [];
N = length(list)
addpath 'C:\Users\feder\Desktop\Tesi\Data\24h Signals'

 %%
for i = N:N
    FOLDER = fullfile(list(i).folder, list(i).name);
    file = dir(FOLDER);
    name = file.name;
    load(name)

    fs_Acc = 64;
    fs_Rot = 64;
    fs_Ecg = 1024; 
    % Filtering - Linear Acceleration (Accelerazione Lineare)
    z_filt = filtHP(Acc.Z_Acc,fs_Acc);
    x_filt = filtroR(Acc.X_Acc,fs_Acc);
    y_filt = filtroR(Acc.Y_Acc,fs_Acc);
    
    
    % Filtering - Angular rate (Velocità Angolare)
    x_rotfilt = filtroR(Rot.X_Rot,fs_Rot);
    y_rotfilt = filtroR(Rot.Y_Rot,fs_Rot);
    z_rotfilt = filtHP(Rot.Z_Rot,fs_Rot);

    Acc_filt = table(x_filt,y_filt,z_filt,Acc.Time_Acc);
    Rot_filt = table(x_rotfilt,y_rotfilt,z_rotfilt,Rot.Time_Rot);

    %save(['FILT-' name ],'Acc_filt','Rot_filt')
    save(['C:\Users\feder\Desktop\Tesi\Data\Filtered Signals\' 'FILT-' name ],'Acc_filt','Rot_filt')
end 

