function ECG_filt = filtECG(ecg,Fs)
%Fs = 1024;                                          % Sampling Frequency
Fn = Fs/2;                                          % Nyquist Frequency
Wp = [1 100]/Fn;                                  % Passband Normalised
Ws = [0.5 105]/Fn;                                  % Stopband Normalised
Rp =  1;                                            % Passband Ripple (Irrelevant in Butterworth)
Rs = 50;                                            % Stopband Attenuation
[n,Wp] = ellipord(Wp,Ws,Rp,Rs);                     % Order Calculation
[z,p,k] = ellip(n,Rp,Rs,Wp);                        % Zero-Pole-Gain 
[sos,g] = zp2sos(z,p,k);                            % Second-Order Section For Stability                               
ECG_filt = filtfilt(sos, g, ecg); 
end 