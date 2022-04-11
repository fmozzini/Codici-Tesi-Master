function  [HR, pos_picchi, amp_picchi,Template] = SCG_template_matching_corr(axis,fs,tempAxis,template)


if isempty(template)
    
    %Zbfilt_REV=axis.^2;
    Zbfilt_REV=axis;


    Zb30sec = {Zbfilt_REV};

    clc
    Template = []; % Matrice contenente gli n Template

    range_pre = 0.20;   % Questi due range rappresentano, moltiplicati per la frequenza, i secondi prima e dopo il picco 

    range_post = 0.20;  % scelto come riferimento.

    % ho solo un blocco di 30 secondi
    blocchi = 1;

    for i=1:blocchi

            Zb10sec_part = Zb30sec{i,1};

            %Zb10sec = Zb30sec_part(1:round(fs*10),1);
            % Zb10sec = Zb30sec_part(round(fs*10):round(fs*20),1);

            [massimo,indice] = max(Zb10sec_part(round(fs*2):round(fs*8),1)); % Finestra da 2 secondi agli 8 secondi. In questo modo si evita la presenza

            indice = indice + round(fs*2)-1;                      % di complessi non completi
            
            Template = cat(1,Template,Zb10sec_part(indice-round(range_pre*fs):indice+round(range_post*fs))); % Template di larghezza 400ms.

    end
else
    Template = template;
    Zbfilt_REV=axis;
    % fede
%      for l = 1:length(Zbfilt_REV)
%             if (Zbfilt_REV(l,1) < 0 )
%                 Zbfilt_REV(l,1) = 0;
%             end 
%      end 
%      %fine fede
    Zb30sec = {Zbfilt_REV};
    blocchi = 1;
    
end

%% Calcolo della Cross-Correlazione
clc


corr_Zbelly = {}; % Matrice contenente le n Cross-Correlazioni



%Calcolo delle n Cross-Correlazioni. La variabile corr_direzione_position

%conterrà i valori di Cross-Correlazione grezzi ossia le Cross-Correlazioni

%senza l'aggiustamento del ritardo. La lunghezza di queste sarà pari alla

%lunghezza dell'intervallo da 30" corrispondente in modo da evitare i

%valori nulli creati dallo zero padding (Il Template è molto più corto del segnale preso

%in considerazione).

%Anche in questo caso il ciclo for esegue N iterazioni pari al numero di

%blocchi



for i=1:blocchi

    Zb30sec_part = Zb30sec{i,1};

    corr_Zbelly_part = xcorr(Template(i,:),Zb30sec_part(:,1));

    corr_Zbelly_part = corr_Zbelly_part(1:length(Zb30sec_part),:);

    corr_Zbelly = cat(1,corr_Zbelly,corr_Zbelly_part);

end

%% Sogliatura della Cross-Correlazione
clc


Soglia = []; % Matrice contentente gli n diversi valori della soglia



% Calcolo della Soglia per ciascun intervallo. Il valore della soglia è

% uguale alla RMS del valore di correlazione di quel dato intervallo 

% aumentato del 50% per evitare i picchi più bassi. Tutti i valori sotto soglia sono posti

% uguali a zero.



for i=1:blocchi

    Soglia = cat(1,Soglia,sqrt(mean((corr_Zbelly{i,1}).^2))*1.5); % Calcolo automatico della soglia (RMS*1.5) 

    for j=1:length(corr_Zbelly{i,1})

        if (corr_Zbelly{i,1}(j,1) < Soglia(i,1))

            corr_Zbelly{i,1}(j,1) = 0;

        end

    end

end


%%
% Rifinitura dei valori massimi della Cross_Correlazione tramite un ciclo

% for che va ad azzerrare tutti i valori non massimi all'interno di un

% piccolo intervallo.
clc


%par_ricerca_cc = 40; % Giustificare la scelta di questo parametro
par_ricerca_cc = 0.4*fs; % SARAH 7.10.2021
% par_ricerca_cc = 0.5*fs; % SARAH 7.10.2021
for r = 1:blocchi

    for i= 1:length(round(corr_Zbelly{r,1}))

        if (corr_Zbelly{r,1}(i,1) > 0) 

            pos = i;

            if (length(corr_Zbelly{r,1}) < pos+par_ricerca_cc) 
                
                [massimo,indice] = max(corr_Zbelly{r,1}(pos:end,1));

                for j= pos:length(corr_Zbelly{r,1})

                    if (corr_Zbelly{r,1}(j,1) ~= massimo)

                        corr_Zbelly{r,1}(j,1)=0;

                    end

                end

                i=i+par_ricerca_cc;

            else

                [massimo,indice] = max(corr_Zbelly{r,1}(pos:pos+par_ricerca_cc,1));

                for j= pos:pos+par_ricerca_cc

                    if (corr_Zbelly{r,1}(j,1) ~= massimo)

                        corr_Zbelly{r,1}(j,1)=0;

                    end

                end

                i=i+par_ricerca_cc;

            end

        end

    end

end

%%
clc
% Eliminazione del primo valore della cross-corr diverso da zero.

for r= 1:blocchi

    for i= 1:length(corr_Zbelly{r,1})

        if (corr_Zbelly{r,1}(i,1) > 0)

            corr_Zbelly{r,1}(i,1) = 0;

            break;

        end

    end

end

% figure
% plot(corr_Zbelly{1, 1}  ) 
% title('Z corr')


%% Unione del segnale e della cross-correlazione

% In questa sezione di codice gli intervalli del segnale così come le cross-correlazioni

% vengono uniti in variabili uniche in modo da evitare la gestione di

% strutture che Zaxispesantiscono il codice.

clc

Zbfilt_REV = []; %Segnale

corr_Zbfilt_REV = [];   %Cross-Correlazione Sogliata



for i= 1:blocchi

    Zbfilt_REV = cat(2,Zbfilt_REV,Zb30sec{i,1});

    corr_Zbfilt_REV = cat(2,corr_Zbfilt_REV,corr_Zbelly{i,1});

end


%% Migliore Correlazione

% Azzeramento di tutti i valori negativi della cross-correlazione se

% presenti
clc


for i=1:length(corr_Zbfilt_REV)

        if(corr_Zbfilt_REV(i,1) <= 0)

            corr_Zbfilt_REV(i,1) = 0;

        end

end


% %%
% % Per evitare di lavorare con tutto il vettore di Cross-Correlazione vengono
% 
% % cancellati tutti i valori uguali a zero in modo da avere solo i valori massimo
% 
% % e una volta fatto questo vengono posti nel vettore punti_iniziali.
% 
% clc
% 
% punti_iniziali = [];
% 
% 
% 
% for i=1:length(corr_Zbfilt_REV)
% 
%     if (corr_Zbfilt_REV(i,1) > 0)
% 
%         punti_iniziali = cat (2,punti_iniziali,i);
% 
%     end
% 
% end
% 
% 
% 
% % Aggiunstamento dei punti iniziali. Per prima cosa viene sottratta la
% 
% % lunghezza del Template e successivamente viene aggiunta 
% 
% % la distanza del picco massimo dal primo punto del Template. In questo
% 
% % modo si ottengono i punti di massima correlazioni relativi al picco
% 
% % voluto.
% 
% 
% 
% punti_iniziali = punti_iniziali - length(Template) + round(range_pre*fs);
% 
% 
% %% Ricerca del primo picco
% 
% % Spesso l'algoritmo non individua correttamente la position del primo
% 
% % picco desiderato quindi si va ad individuarlo in maniera manuale.
% 
% clc
% 
% [massimo,indice] = max(Zbfilt_REV(1:punti_iniziali(1)-20,1));
% 
% punti_iniziali = [indice, punti_iniziali];
% 
% 
% 
% if (punti_iniziali(1,end) > length(Zbfilt_REV))
% 
%     punti_iniziali(1,end) = length(Zbfilt_REV);
% 
% end
% 
% 
% 
% %% Ricerca dei massimi a partire dai valori di cross-correlazione
% 
% % In questa sezione vengono ricercate in prima Zaxisprossimazione le
% 
% % posizioni delle waves desiderate. Si cerca un valore massimo all'interno di un certo
% 
% % range a partire dai punti iniziali trovati.
% 
% clc
% 
% J_prima_ricerca = []; 
% 
% 
% 
% par_prima_ricerca = round(mean(diff(punti_iniziali))/2);
% 
% 
% 
% % La ricerca del picco massimo viene fatta in un intervallo di n millisecondi.
% 
% % La scelta del range viene fatta in modo da avere un intervallo totale pari 
% 
% % alla media della durata del complesso. Questo intervallo è soggetto dipendente. 
% 
% % I passi successivi dell'algoritmo aggiusterranno i possibili errori generati. 
% 
% 
% 
% for i=1:length(punti_iniziali)
% 
%     if(punti_iniziali(i)+par_prima_ricerca > length(Zbfilt_REV))
% 
%         break;
% 
%     elseif(punti_iniziali(i)-par_prima_ricerca < 0)
% 
%         [massimo,indice] = max(Zbfilt_REV(1:punti_iniziali(i)+par_prima_ricerca,1));
% 
%         J_prima_ricerca = cat(2,J_prima_ricerca,indice);
% 
%     else[massimo,indice] = max(Zbfilt_REV(punti_iniziali(i)-par_prima_ricerca:punti_iniziali(i)+par_prima_ricerca,1)); 
% 
%         J_prima_ricerca = cat(2,J_prima_ricerca,indice+punti_iniziali(i)-par_prima_ricerca-1);
% 
%     end
% 
% end
% 
% 
% 
% % Il ciclo for e la successiva riga di codice eliminano i doppioni di punti di waves trovati
% 
% for i=1:length(J_prima_ricerca)-1
% 
%     if(J_prima_ricerca(i+1)==J_prima_ricerca(i))
% 
%         J_prima_ricerca(i)=0;
% 
%     end
% 
% end
% 
% 
% 
% J_prima_ricerca(J_prima_ricerca==0) = [];
% 
% %% Ricerca waves troppo lontane
% 
% % Una volta trovate le posizioni delle waves si procede con una lunga fase
% 
% % di correzione della position di queste.
% 
% % Per prima cosa si vanno a cercare tutte quelle waves che si
% 
% % trovano ad una distanza superiore al 50% rispetto alla distanza media.
% 
% % Questa distanza viene aggiornata ad ogni ciclo for in modo da avere
% 
% % sempre un valore attuale che rende la correzione più precisa.
% 
% clc
% 
% distanza_media = round(mean(diff(J_prima_ricerca))); % Calcolo della distanza media tra i picchi J 
% 
% 
% 
% par_ricerca_waves_mancanti = round(distanza_media/3);
% 
% % La ricerca del picco mancante viene fatta nell'intervallo compreso tra i punti ritrovati
% 
% % iniziando però da  distanza media picchi/3 campioni dopo la prima onda e dello stesso ammontare prima. In
% 
% % questo modo si evita di trovare un'onda Zaxispartenente ai complessi considerati
% 
% 
% 
% J_mancanti = [];
% 
% J_seconda_ricerca = J_prima_ricerca;
% 
% 
% 
% for i= 1:length(J_prima_ricerca)-1
% 
%     if (J_prima_ricerca(i+1)-J_prima_ricerca(i) > round(distanza_media*1.5))
% 
%         [massimo,indice] = max(Zbfilt_REV(J_prima_ricerca(1,i)+par_ricerca_waves_mancanti:J_prima_ricerca(1,i+1)-par_ricerca_waves_mancanti,1));
% 
%         J_seconda_ricerca = cat(2,J_seconda_ricerca,indice+J_prima_ricerca(1,i)+par_ricerca_waves_mancanti-1);
% 
%         J_seconda_ricerca = sort(J_seconda_ricerca);
% 
%         distanza_media = round(mean(diff(J_seconda_ricerca)));
% 
%     end
% 
% end
% 
% 
% 
% % Fatta una prima ricerca delle waves mancanti il processo viene ripeturo i
% 
% % in modo da correggere eventuali waves non trovate. In questo caso
% 
% % la media non viene aggiornata poichè i cambiamenti sarebbero lievi.
% 
% 
% 
% for i= 1:length(J_seconda_ricerca)-1
% 
%     if (J_seconda_ricerca(i+1)-J_seconda_ricerca(i) > round(distanza_media*1.5))
% 
%         [massimo,indice] = max(Zbfilt_REV(J_seconda_ricerca(1,i)+par_ricerca_waves_mancanti:J_seconda_ricerca(1,i+1)-par_ricerca_waves_mancanti,1));
% 
%         J_mancanti = cat(2,J_mancanti,indice+J_seconda_ricerca(1,i)+par_ricerca_waves_mancanti-1);
% 
%     end
% 
% end
% 
% 
% 
% % Le waves trovate vengono unite e ordinate.
% 
% 
% 
% J_seconda_ricerca = cat(2,J_seconda_ricerca,J_mancanti);
% 
% J_seconda_ricerca = sort(J_seconda_ricerca);
% 
% 
% 
% 
% %% Ricerca waves troppo ravvicinate
% 
% % In questa sezione si lavora similmente a prima solamente che vengono
% 
% % ricercate le waves che sono ad una distanza minore del 50% rispetto la
% 
% % media
% 
% clc
% 
% distanza_media = round(mean(diff(J_seconda_ricerca))); 
% 
% Aggiungi_J = [];
% 
% 
% 
% % Il ciclo for che svolge questa funzione viene fatto per due volte in modo
% 
% % da correggere eventuali errori.
% 
% 
% 
% for cicle = 1:2
% 
%     for i=1:length(J_seconda_ricerca)-1
% 
%         if (J_seconda_ricerca(i+1)-J_seconda_ricerca(i) < distanza_media*0.5)
% 
%             [massimo,indice] = max(Zbfilt_REV(J_seconda_ricerca(1,i):J_seconda_ricerca(1,i+1),1));
% 
%             Aggiungi_J = cat(2,Aggiungi_J,indice+J_seconda_ricerca(1,i)-1);
% 
%         end
% 
%     end
% 
% end
% 
% 
% 
% % Una volta fatto questo viene creata una variabile waves_J che contiene i
% 
% % valori delle waves J indentificate in POS2 Zaxis.
% 
% 
% 
% waves_J = [];
% 
% j = 1;
% 
% for i= 1:length(J_seconda_ricerca)
% 
%     if j > length(Aggiungi_J)
% 
%         break;
% 
%     elseif J_seconda_ricerca(1,i) < Aggiungi_J(1,j) && J_seconda_ricerca(1,i+1) > Aggiungi_J (1,j)
% 
%         waves_J = cat(2,waves_J, Aggiungi_J(1,j));
% 
%         i = i+1;
% 
%         j = j+1;
% 
%     else
% 
%         waves_J = cat(2,waves_J, J_seconda_ricerca(1,i));
% 
%     end
% 
% end
% 
% 
% 
% waves_J = cat(2,waves_J, J_seconda_ricerca(1,i:end));
% 
% 
% 
% j=1;
% 
% for i= 1:length(waves_J)
% 
%     if j > length(Aggiungi_J)
% 
%     elseif (waves_J(1,i) == Aggiungi_J (1,j))
% 
%         waves_J(1,i+1) = 0;
% 
%         j = j+1;
% 
%     end
% 
% end
% 
% 
% 
% waves_J(waves_J==0) = [];
% 
% 
% 
% for i=1:length(waves_J)-1
% 
%     if(waves_J(i+1)==waves_J(i))
% 
%         waves_J(i)=0;
% 
%     end
% 
% end
% 
% 
% 
% waves_J(waves_J==0) = [];
% 
% 
% 
% %% Correzione dalla Serie
% 
% % Una volta eseguiti gli aggiustamenti viene calcolata la media della
% 
% % serie_JJ e la serie parziale.
% 
% clc
% 
% media_serie_JJ = mean(diff(waves_J));
% 
% serie_parziale = diff(waves_J);
% 
% 
% 
% % Vengono effettuate correzioni similmente a prima sui valori delle serie
% 
% % Funzionano allo stesso modo delle precedenti però evitano errori
% 
% 
% 
% par_serie = round(media_serie_JJ/3);
% 
% 
% 
% % Si parte con la correzione su tempi molto lunghi
% 
% for i=1:length(1:length(serie_parziale))
% 
%        if (serie_parziale(i) > media_serie_JJ*1.5)
% 
%            [massimo,indice] = max(Zbfilt_REV(waves_J(i)+par_serie:waves_J(i+1)-par_serie,1));
% 
%            waves_J = cat(2,waves_J,indice+waves_J(i)+par_serie-1);
% 
%        end
% 
% end
% 
% 
% 
% waves_J = sort(waves_J);
% 
% 
% 
% media_serie_JJ = mean(diff(waves_J));
% 
% serie_parziale = diff(waves_J);
% 
% 
% 
% % Si effettua la stessa cosa per durate molto brevi
% 
% for i=1:length(waves_J)-1
% 
%     if (serie_parziale(i) < media_serie_JJ*0.5)
% 
%         [massimo,indice] = max(Zbfilt_REV(waves_J(i):waves_J(i+1),1));
% 
%         indice = indice + waves_J(i)-1;
% 
%         waves_J(i) = 0;
% 
%         waves_J(i+1)= indice;
% 
%     end
% 
% end
% 
% 
% 
% waves_J(waves_J==0) = [];
% 
% waves_J = sort(waves_J);
% 
% 
% 
% 
% %% Correzione waves
% 
% % Per la direzione Zaxis in position 2 (Belly) può cZaxisitare che l'algoritmo
% 
% % individui l'onda L invece che la J e quindi è opportuno effettuare una
% 
% % correzione che indivudui l'onda J. Per fare questo è possibile cercare in
% 
% % un intervallo a sinistra dell'onda trovata un'altra onda di dimensione
% 
% % comparabile a quella trovata. In questo modo se viene trovata l'onda L si
% 
% % ritorna alla J mentre se viene individuata l'onda J si rimane in questo
% 
% % punto poichè l'onda H non ha mai misura comparabile alla J.Il criterio di
% 
% % comparazione delle waves è maggiore o almeno il 70%.
% 
% clc
% 
% 
% par_Jbasse = round(mean(diff(waves_J))/5);
% 
% 
% 
% for i=1:length(waves_J)
% 
%     if (waves_J(1,i)-par_Jbasse <= 0)
% 
%         part = Zbfilt_REV(1:waves_J(1,i)-5,1);
% 
%         [massimo,indice] = max(part);
% 
%         valore_comp = Zbfilt_REV(waves_J(1,i),1);
% 
%         indice = indice + waves_J(1,i)-1;
% 
%         if (massimo > valore_comp*0.7)
% 
%             waves_J (1,i) = indice;
% 
%         end
% 
%     else
% 
%         part = Zbfilt_REV(waves_J(1,i)-par_Jbasse:waves_J(1,i)-5,1);
% 
%         [massimo,indice] = max(part);
% 
%         valore_comp = Zbfilt_REV(waves_J(1,i),1);
% 
%         indice = indice + waves_J(1,i)-par_Jbasse-1;
% 
%         if (massimo > valore_comp*0.7)
% 
%             waves_J (1,i) = indice;
% 
%         end
% 
%     end
% 
% end
% 
% 
% 
% waves_J_Zaxis_REV = waves_J;
% 
% 















































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
close all

% figure(1)
% plot (Template)
% savefig(strcat('Template',tempAxis,'.fig'))
% close all
%%

corr_Zbfilt_dritto = corr_Zbfilt_REV(end:-1:1,1);
picco = find(corr_Zbfilt_dritto~=0);
picco

% soglia tempo
picco_temp = picco;
count = 0;

for i = 1:size(picco_temp)-1
    if strcmp (picco_temp(i), NaN) ==1
        continue
    end
% CAMBIAMENTO FEDE 
    if (picco_temp(i+1,1) - picco_temp(i,1)) < 0.4*fs 
%     if (picco_temp(i+1,1) - picco_temp(i,1)) < 0.5*fs
        picco_temp (i+1) = NaN;
    end
end
picco_temp = picco_temp(isnan(picco_temp)==0);

% figure
% plot (axis,'k')
% hold on
% plot (picco_temp, axis(picco_temp),'r*')


% 
% figure
% plot (axis,'k')
% hold on
% plot (picco_temp, axis(picco_temp),'r*')
% savefig(strcat('Picchi_TempMatch',tempAxis,'.fig'))
% % pause CON PAUSE OGNI VOLTA DEVO SCHIACCHIARE ENTER 
% close all
% 
pos_picchi = picco_temp;
amp_picchi = axis(picco_temp);
pos_picchi(isnan(amp_picchi)) = [];
amp_picchi(isnan(amp_picchi)) = [];


%% Compute HR
clc
HR = 60./(diff(pos_picchi))*fs;

%%




















end

