%% Kinetic Energy 
clear all
close all
clc

%% 
folder = 'C:\Users\feder\Desktop\Tesi\Data\Filtered Signals';
list = dir(folder);
list(1) = [];
list(1) = [];
list(21) = [];
list(21) = [];
N = length(list)

data = xlsread('C:\Users\feder\Desktop\Tesi\Info Pazienti.xlsx','E:G');

%%
for i = 1:N
    FOLDER = fullfile(list(i).folder, list(i).name);
    file = dir(FOLDER);
    name = file.name;
    load(name)
    
    weight = data(i,2)
    heigth = data(i,3)
    %% Kinetic Energy 
    % Lineare
    cinetical = KE(weight,Acc_filt.x_filt,Acc_filt.y_filt,Acc_filt.z_filt,Acc_filt.Var4);
    
    % Rotazionale 
    cineticar = KEr(weight,heigth,Rot_filt.x_rotfilt,Rot_filt.y_rotfilt,Rot_filt.z_rotfilt);
   
    name_KE = erase(name,"FILT-")
    save(['C:\Users\feder\Desktop\Tesi\Data\LinearKE\' 'LinearKE-' name_KE],'cinetical')
    save(['C:\Users\feder\Desktop\Tesi\Data\RotationalKE\' 'RotationalKE-' name_KE],'cineticar')
    %save(['DatiKE-' name_KE],'weight','heigth')
    
end 

