% View data
clear all
close all
clc
%%
load '2022-01-18 20.30.02 Subject 03.mat'             
% Lorenzo Monti -> Peso: 83 Kg, Altezza: 188 cm
% Inizio: 20:00, Sonno: 1:30-8:30
% 0: 20-21; 1: 21-22; 2: 22-23; 3: 23-24; 4: 24-1;
% 5: 1-2; 6: 2-3; 7:3-4; 8: 4-5; 9: 5-6; 10: 6-7;
% 11: 7-8; 12: 8-9
% Dorme da 5 ora a 12

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


a1 = nexttile
plot(x_filt/64)
a2 = nexttile
plot(x_rotfilt/64)
a3 = nexttile
plot(Ecg.Values/1024)
linkaxes([a1 a2 a3],'x')

%% 10 secondi
Hour = 3600*64;
Ecg_Hour = 3600*1024;
% Guardo 10 secondi durante il sonno

figure()
subplot(211)
plot(Acc.X_Acc(6*Hour:6*Hour+640),'r')
hold on;
plot(x_filt(6*Hour:6*Hour+640),'b')
subplot(212)
plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+10240))
sgtitle('10 secondi notte (2:00) Acc X ')

figure()
subplot(211)
plot(Acc.Y_Acc(6*Hour:6*Hour+640),'r')
hold on;
plot(y_filt(6*Hour:6*Hour+640),'b')
subplot(212)
plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+10240))
sgtitle('10 secondi notte (2:00) Acc Y ')

figure()
subplot(211)
plot(Acc.Z_Acc(6*Hour:6*Hour+640),'r')
hold on;
plot(z_filt(6*Hour:6*Hour+640),'b')
subplot(212)
plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+10240))
sgtitle('10 secondi notte (2:00) Acc Z ')

figure()
subplot(211)
plot(Rot.X_Rot(6*Hour:6*Hour+640),'r')
hold on;
plot(x_rotfilt(6*Hour:6*Hour+640),'b')
subplot(212)
plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+10240))
sgtitle('10 secondi notte (2:00) Rot X ')

figure()
subplot(211)
plot(Rot.Y_Rot(6*Hour:6*Hour+640),'r')
hold on;
plot(y_rotfilt(6*Hour:6*Hour+640),'b')
subplot(212)
plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+10240))
sgtitle('10 secondi notte (2:00) Rot Y ')

figure()
subplot(211)
plot(Rot.Z_Rot(6*Hour:6*Hour+640),'r')
hold on;
plot(z_rotfilt(6*Hour:6*Hour+640),'b')
subplot(212)
plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+10240))
sgtitle('10 secondi notte (2:00) Rot Z ')

%% 5 secondi
Hour = 3600*64;
Ecg_Hour = 3600*1024;
% Guardo 10 secondi durante il sonno

figure()
subplot(211)
plot(Acc.X_Acc(10*Hour:10*Hour+320),'r')
hold on;
plot(x_filt(10*Hour:10*Hour+320),'b')
subplot(212)
plot(Ecg.Values(10*Ecg_Hour:10*Ecg_Hour+5120))
sgtitle('5 secondi notte (2:00) Acc X ')

figure()
subplot(211)
plot(Acc.Y_Acc(6*Hour:6*Hour+320),'r')
hold on;
plot(y_filt(6*Hour:6*Hour+320),'b')
subplot(212)
plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+5120))
sgtitle('5 secondi notte (2:00) Acc Y ')

figure()
subplot(211)
plot(Acc.Z_Acc(6*Hour:6*Hour+320),'r')
hold on;
plot(z_filt(6*Hour:6*Hour+320),'b')
subplot(212)
plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+5120))
sgtitle('5 secondi notte (2:00) Acc Z ')

figure()
subplot(211)
plot(Rot.X_Rot(10*Hour:10*Hour+320),'r')
hold on;
plot(x_rotfilt(10*Hour:10*Hour+320),'b')
subplot(212)
plot(Ecg.Values(10*Ecg_Hour:10*Ecg_Hour+5120))
sgtitle('5 secondi notte (2:00) Rot X ')

figure()
subplot(211)
plot(Rot.Y_Rot(6*Hour:6*Hour+320),'r')
hold on;
plot(y_rotfilt(6*Hour:6*Hour+320),'b')
subplot(212)
plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+5120))
sgtitle('5 secondi notte (2:00) Rot Y ')

figure()
subplot(211)
plot(Rot.Z_Rot(6*Hour:6*Hour+320),'r')
hold on;
plot(z_rotfilt(6*Hour:6*Hour+320),'b')
subplot(212)
plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+5120))
sgtitle('5 secondi notte (2:00) Rot Z ')
%% Kinetic Energy 
% Lineare
cinetical = KE(53,x_filt,y_filt,z_filt,Acc.Time_Acc);

figure()
subplot(211)
plot(cinetical(10*Hour:10*Hour+630)),title('Linear Kinetic Energy 2:00 per 10 sec')
subplot(212)
plot(Ecg.Values(10*Ecg_Hour:10*Ecg_Hour+10240)),title('Ecg 2:00 per 10 sec')

%%
% Rotazionale
cineticar = KEr(53,165,x_rotfilt,y_rotfilt,z_rotfilt);

figure()
subplot(211)
plot(cineticar(10*Hour:10*Hour+630)),title('Rotational Kinetic Energy 2:00 per 10 sec')
subplot(212)
plot(Ecg.Values(10*Ecg_Hour:10*Ecg_Hour+10240)),title('Ecg 2:00 per 10 sec')



%% PAN TOMPKIN
Hour = 3600*64;
Ecg_Hour = 3600*1024;
[qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+5120),1024,1)

























