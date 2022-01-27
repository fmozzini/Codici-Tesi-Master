% Read files 
clear all
close all
clc
 
clear filename
filename = strcat('/Users/feder/Desktop/Tesi/Movisens acquisitions/2021-06-12 09.38.41 MAIOLO/',...
    'acc.csv'); % INSERIRE PATH CORRETTO
acc = readtable(filename);
acc.Properties.VariableNames = {'x' 'y' 'z'};
fs_acc = 64; % Hz
 
clear filename
filename = strcat('/Users/feder/Desktop/Tesi/Movisens acquisitions/2021-06-12 09.38.41 MAIOLO/',...
    'angularrate.csv'); % INSERIRE PATH CORRETTO
gyro = readtable(filename);
gyro.Properties.VariableNames = {'x' 'y' 'z'};
fs_gyro = 64; % Hz
 
clear filename
filename = strcat('/Users/feder/Desktop/Tesi/Movisens acquisitions/2021-06-12 09.38.41 MAIOLO/',...
    'ecg.csv'); % INSERIRE PATH CORRETTO
ecg = readtable(filename);
ecg.Properties.VariableNames = {'signal'};
clear filename
fs_ecg = 1024; % Hz

time_Acc = 1:size(acc,1); %5529584 -> ho un valore ogni 0.0156 secondi
% ogni secondo ho 64 valori -> 24 ore 
% 86556 sec/60 min/60 -> ore
time_Ecg = 1:size(ecg,1); %88473343 -> ho un valore ogni 1/1024 secondi
% ogni secondo ho 1024 valori -> 24 ore


% SBAGLIATO, NON SONO SECONDI MA HZ - campioni
% plot Acc x,y,z
timeA_fs = time_Acc./fs_acc;
f1 = figure()
subplot(311)
sgtitle('Acceleration x,y,z')
plot(time_Acc./fs_acc, acc.x,'r'), xlabel('Time [s]'), ylabel('[ m/s^2 ]'),legend('acc(x)'),axis tight
subplot(312)
plot(time_Acc./fs_acc, acc.y,'g'), xlabel('Time [s]'), ylabel('[ m/s^2 ]'),legend('acc(y)'),axis tight
subplot(313)
plot(time_Acc/fs_acc, acc.z,'b'), xlabel('Time [s]'), ylabel('[ m/s^2 ]'),legend('acc(z)'),axis tight

% plot Gyro x,y,z
f2 = figure()
subplot(311)
sgtitle('Angular rate x,y,z')
plot(time_Acc, gyro.x,'r'), xlabel('Time [s]'), ylabel('[ rad/s ]'),legend('gyro(x)'),axis tight
subplot(312)
plot(time_Acc, gyro.y,'g'), xlabel('Time [s]'), ylabel('[ rad/s ]'),legend('gyro(y)'),axis tight
subplot(313)
plot(time_Acc, gyro.z,'b'), xlabel('Time [s]'), ylabel('[ rad/s ]'),legend('gyro(z)'),axis tight

% plot ECG
f3 = figure()
plot(time_Ecg, ecg.x3747), xlabel('Time [s]'), legend('ecg'), axis tight
title('Ecg')

for i = 0:10
    figure()
    plot(time_Acc((640*i+1):(640*(i+1))),acc.x((640*i+1):(640*(i+1))))
end 

time_Acc = (0:numel(acc.x)-1)/fs_acc;
time_AccH = time_Acc/3600;
time_Gy = (0:numel(gyro.x)-1)/fs_gyro;
time_GyH = time_Gy/3600;
time_Ecg = (0:numel(ecg.x3747)-1)/fs_ecg;
time_EcgH = time_Ecg/3600;
time_AccM = time_Acc/60;
time_GyM = time_Gy/60;
time_EcgM = time_Ecg/60;

figure
ax(1) = subplot(3,1,1);
plot(time_AccH,acc.x,'k')
ylabel('Accelerometer x')
grid on
ax(2) = subplot(3,1,2); 
plot(time_GyH,gyro.x,'r')
ylabel('Gyroscope x')
grid on
ax(3) = subplot(3,1,3); 
plot(time_EcgH,ecg.x3747)
ylabel('ECG')
grid on
xlabel('Time (Hours)')
linkaxes(ax(1:3),'x')
axis([0 24 -1 1])
 
for i = 2:23
    figure()
    plot(time_AccH(230400*(i-1):230400*i),acc.x(230400*(i-1):230400*i),'k')
end

% quanti valori ho in 1 ora? ho 1024 Hz 
% 1024 valori al secondo -> x 60 x 60 -> 1:3686400 (quanti all'ora)
plot(time_EcgH(1:3686400),ecg.x3747(1:3686400))
ylabel('ECG')
%axis([0 24 250 4500])

% 1024 valori al secondo -> x 60 = 61440
plot(time_Ecg(1:61440),ecg.x3747(1:61440))
ylabel('ECG')

% 1024 al secondo -> 10 secondi
f1 = figure()
plot(time_Ecg(1:10*1024),ecg.x3747(1:10*1024))
ylabel('ECG')
f2 = figure()
plot(time_Ecg(1:20*1024),ecg.x3747(1:20*1024))
f3 = figure()
plot(time_Ecg(10*1024:20*1024),ecg.x3747(10*1024:20*1024))

% Analisi su finestra di 10 sec
figure
ax(1) = subplot(3,1,1);
plot(time_Acc(1:10*64),acc.x(1:10*64),'k')
ylabel('Accelerometer x')
grid on
ax(2) = subplot(3,1,2); 
plot(time_Gy(1:10*64),gyro.x(1:10*64),'r')
ylabel('Gyroscope x')
grid on
ax(3) = subplot(3,1,3); 
plot(time_Ecg(1:10*1024),ecg.x3747(1:10*1024))
ylabel('ECG')
grid on
xlabel('Time (secs)')
linkaxes(ax(1:3),'x')
%axis([0 24 -1 1])

% Analisi su 1 sec
for i =  0:10
figure()
ax(1) = subplot(3,1,1);
plot(time_Acc((64*10*i+1):(64*10*(i+1))),acc.z((64*10*i+1):(64*10*(i+1))),'k')
ylabel('Accelerometer z')
axis([10*i+1 10*(i+1) -2050 3450])
grid on 
ax(2) = subplot(3,1,2); 
plot(time_Gy((64*10*i+1):(64*10*(i+1))),gyro.z((64*10*i+1):(64*10*(i+1))),'r')
ylabel('Gyroscope z')
axis([10*i+1 10*(i+1) -1670 1550])
grid on
ax(3) = subplot(3,1,3); 
plot(time_Ecg((1024*10*i+1):(1024*10*(i+1))),ecg.x3747((1024*10*i+1):(1024*10*(i+1))))
ylabel('ECG')
axis([10*i+1 10*(i+1) 0 4100])
grid on
xlabel('Time (secs)')
linkaxes(ax(1:3),'x')
end 

figure()
plot(time_Ecg(1:1024*20),ecg.x3747(1:1024*20))