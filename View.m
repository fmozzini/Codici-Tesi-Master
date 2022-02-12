% View data
clear all
close all
clc
%%
load '2022-02-09 20.00.01 Subject 09'
%% Check dimension
Ore_SCG = (size(Acc.X_Acc,1))/(64*3600)
Ore_ECG = (size(Ecg.Values,1))/(1024*3600)
%% Check PT 
load 'PT-2021-01-16 15.01.53.mat'
% figure()
% plot(ecg_FILT)
% hold on 
% scatter(qrs_I,qrs_AMP,'m')

wind = 10*1024;
for i = 1:10
    if i == 1
        figure()
        plot((i-1)*wind+1:i*wind,ecg_FILT((i-1)*wind+1:i*wind))
        hold on
        scatter(qrs_I(peaksFORwindow(i):peaksFORwindow(i)),qrs_AMP(peaksFORwindow(i):peaksFORwindow(i)),'m')
    else 
        figure()
        plot((i-1)*wind+1:i*wind,ecg_FILT((i-1)*wind+1:i*wind))
        hold on
        scatter(qrs_I(peaksFORwindow(i-1)+1:peaksFORwindow(i)),qrs_AMP(peaksFORwindow(i-1)+1:peaksFORwindow(i)),'m')
    end 
end 
%%
% Quanti battiti ci sono al minuto (60 sec) e quanto in ogni finestra di 10 s?
n_finestre = size(peaksFORwindow,1);
n_picchi_min = zeros(round(n_finestre/6),1);
n_picchi_10s = zeros(n_finestre,1);
for i = 1:n_finestre-1 
    if i == 1
        n_picchi_10s(i) = peaksFORwindow(i);
    else 
        n_picchi_10s(i) = peaksFORwindow(i)-peaksFORwindow(i-1);
    end
end 

mean_picchi_10s = mean(n_picchi_10s)
mean_picchi_min = 6*mean_picchi_10s    

% %% 10 secondi
% Hour = 3600*64;
% Ecg_Hour = 3600*1024;
% % Guardo 10 secondi durante il sonno
% 
% figure()
% subplot(211)
% plot(Acc.X_Acc(6*Hour:6*Hour+640),'r')
% hold on;
% plot(x_filt(6*Hour:6*Hour+640),'b')
% subplot(212)
% plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+10240))
% sgtitle('10 secondi notte (2:00) Acc X ')
% 
% figure()
% subplot(211)
% plot(Acc.Y_Acc(6*Hour:6*Hour+640),'r')
% hold on;
% plot(y_filt(6*Hour:6*Hour+640),'b')
% subplot(212)
% plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+10240))
% sgtitle('10 secondi notte (2:00) Acc Y ')
% 
% figure()
% subplot(211)
% plot(Acc.Z_Acc(6*Hour:6*Hour+640),'r')
% hold on;
% plot(z_filt(6*Hour:6*Hour+640),'b')
% subplot(212)
% plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+10240))
% sgtitle('10 secondi notte (2:00) Acc Z ')
% 
% figure()
% subplot(211)
% plot(Rot.X_Rot(6*Hour:6*Hour+640),'r')
% hold on;
% plot(x_rotfilt(6*Hour:6*Hour+640),'b')
% subplot(212)
% plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+10240))
% sgtitle('10 secondi notte (2:00) Rot X ')
% 
% figure()
% subplot(211)
% plot(Rot.Y_Rot(6*Hour:6*Hour+640),'r')
% hold on;
% plot(y_rotfilt(6*Hour:6*Hour+640),'b')
% subplot(212)
% plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+10240))
% sgtitle('10 secondi notte (2:00) Rot Y ')
% 
% figure()
% subplot(211)
% plot(Rot.Z_Rot(6*Hour:6*Hour+640),'r')
% hold on;
% plot(z_rotfilt(6*Hour:6*Hour+640),'b')
% subplot(212)
% plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+10240))
% sgtitle('10 secondi notte (2:00) Rot Z ')
% 
% %% 5 secondi
% Hour = 3600*64;
% Ecg_Hour = 3600*1024;
% % Guardo 10 secondi durante il sonno
% 
% figure()
% subplot(211)
% plot(Acc.X_Acc(10*Hour:10*Hour+320),'r')
% hold on;
% plot(x_filt(10*Hour:10*Hour+320),'b')
% subplot(212)
% plot(Ecg.Values(10*Ecg_Hour:10*Ecg_Hour+5120))
% sgtitle('5 secondi notte (2:00) Acc X ')
% 
% figure()
% subplot(211)
% plot(Acc.Y_Acc(6*Hour:6*Hour+320),'r')
% hold on;
% plot(y_filt(6*Hour:6*Hour+320),'b')
% subplot(212)
% plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+5120))
% sgtitle('5 secondi notte (2:00) Acc Y ')
% 
% figure()
% subplot(211)
% plot(Acc.Z_Acc(6*Hour:6*Hour+320),'r')
% hold on;
% plot(z_filt(6*Hour:6*Hour+320),'b')
% subplot(212)
% plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+5120))
% sgtitle('5 secondi notte (2:00) Acc Z ')
% 
% figure()
% subplot(211)
% plot(Rot.X_Rot(10*Hour:10*Hour+320),'r')
% hold on;
% plot(x_rotfilt(10*Hour:10*Hour+320),'b')
% subplot(212)
% plot(Ecg.Values(10*Ecg_Hour:10*Ecg_Hour+5120))
% sgtitle('5 secondi notte (2:00) Rot X ')
% 
% figure()
% subplot(211)
% plot(Rot.Y_Rot(6*Hour:6*Hour+320),'r')
% hold on;
% plot(y_rotfilt(6*Hour:6*Hour+320),'b')
% subplot(212)
% plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+5120))
% sgtitle('5 secondi notte (2:00) Rot Y ')
% 
% figure()
% subplot(211)
% plot(Rot.Z_Rot(6*Hour:6*Hour+320),'r')
% hold on;
% plot(z_rotfilt(6*Hour:6*Hour+320),'b')
% subplot(212)
% plot(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+5120))
% sgtitle('5 secondi notte (2:00) Rot Z ')
% %% Kinetic Energy 
% % Lineare
% cinetical = KE(53,x_filt,y_filt,z_filt,Acc.Time_Acc);
% 
% figure()
% subplot(211)
% plot(cinetical(10*Hour:10*Hour+630)),title('Linear Kinetic Energy 2:00 per 10 sec')
% subplot(212)
% plot(Ecg.Values(10*Ecg_Hour:10*Ecg_Hour+10240)),title('Ecg 2:00 per 10 sec')
% 
% %%
% % Rotazionale
% cineticar = KEr(53,165,x_rotfilt,y_rotfilt,z_rotfilt);
% 
% figure()
% subplot(211)
% plot(cineticar(10*Hour:10*Hour+630)),title('Rotational Kinetic Energy 2:00 per 10 sec')
% subplot(212)
% plot(Ecg.Values(10*Ecg_Hour:10*Ecg_Hour+10240)),title('Ecg 2:00 per 10 sec')
% 
% 
% 
% %% PAN TOMPKIN
% Hour = 3600*64;
% Ecg_Hour = 3600*1024;
% [qrs_amp_raw,qrs_i_raw,delay]=pan_tompkin(Ecg.Values(6*Ecg_Hour:6*Ecg_Hour+5120),1024,1)

























