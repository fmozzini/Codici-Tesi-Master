function sigfilt=filtroR(signal,fs)

[b,a]=butter(4,1/(fs/2),'high');
%    [b,a]=butter(8,0.02,'high');
 sigfilt=filtfilt(b,a,signal);
[b,a]=butter(4,30/(fs/2),'low');
%    [b,a]=butter(8,0.6,'low');
 sigfilt=filtfilt(b,a,sigfilt);
 
end



