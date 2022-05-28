function [HR_min,HR_5min] = HR(qrs_I,fs_ECG)
    qrs_min =  (qrs_I./fs_ECG)./60; 
    n_cinqueminuti = round(qrs_min(end)/5);
    n_min = round(qrs_min(end));
    HR_5min = zeros(n_cinqueminuti,1);
    HR_min = zeros(n_min,1);
    for i = 1:n_min
        if i ==1
            posizione = find(qrs_min<i);
            valori = qrs_min(posizione);
        else
            posizione1 = find(qrs_min<(i+1));
            valori1 = qrs_min(posizione1);
            posizione = find(valori1>i);
            valori = valori1(posizione);
        end
        HR_min(i) = length(valori);
    end
    m = 0;
    for i = 1:5:n_min-5
        m = m+1;
        HR_5min(m) = mean(HR_min(i:i+4));
    end 
 end