function tw(action)
% Time Warping Analysis

% ver 1.8 - 13/3/2014
% modified loading to select BCG channels from multichannel recordings
% ver 1.7 - 6/7/2001
% incluso 'TWMIO1NN' per lettura file originale e creazione matrici con waveforms

% ver 1.6 - 23/5/00
% aggiunto calcolo IVRT (punti iniziali e finali, durata) con correzione/stampa(.xls .mat)
% ver 1.5 - 19/5/00
% aggiunto calcolo e stampa in .xls dei parametri atrfil,tPER e tPAFR
% corretto calcolo time13 e amp13
% corretto opzione correzione pafr
% introdotto opzione correzione pfr e per
% ver 1.4 - 8/5/00
% calcolo parametri LAA
% ver 1.3 - 28/3/00
% corretto estrazione punti fiduciari e salvataggio valori estratti
% aggiunto opzione 'computeTW'(calcolo matrici WF) e 'serie'(estrazione serie di variabilità)
% DA VEDERE !!!
% ver 1.2 - 15/3/00
% aggiunto opzione 'loadAW' e 'extract' con salvataggio dei parametri su foglio Excel
% ver 1.1 - 14/3/00
% aggiunto opzione 'ave' e sistemato il 'load' con selezione battiti iniziale e finale
% ver. 1.0 - 21/02/2000

global NOME LAA matorig mat1 n M mat2 nWF NWF nAW NAW Temp carico caricoAW  caricoPAR caricoWF T template WFX WFY tedv tkneeD2 tpfr tknee2D2 tpafr tper Fc PERC111
% matorig: matrice originale con tutti i battiti
% mat1	: matrice con solo i cicli relativi al calcolo del AW
% n		: # righe di mat1
% M		: # colonne di mat1
% mat2	: matrice con solo i cicli relativi al calcolo delle WF
% nWF		: # battito iniziale rispetto a matorig per il calcolo delle WF
% NWF		: # battito finale rispetto a matorig per il calcolo delle WF
% Temp	: matrice contenente le varie AW estratte
% carico		: 1 se ho letto il file con i vari cicli
% caricoAW	: 1 se ho calcolato il template
% caricoPAR	: 1 se ho calcolato i parametri del template
% caricoWF	: 1 se ho calcolato le matrici delle WF rispetto al template
% T		: lunghezza del template
% template: AW calcolata
% WFX e WFY: matrici delle WF calcolate con TWDYN0
Fc=100;
LAA=0;

switch (action)
       
    case 'load',
        [NOME,PERC1]=uigetfile('*.mat','Import File with waveforms to be analyzed');
        if NOME==0,
            warndlg('No File selected !','!! WARNING !!');
        else
            ST0=['PERC','FILE'];
            STR=strrep(ST0,'PERC',PERC1);
            STR1=strrep(STR,'FILE',NOME);
            % Leggo file
            load(STR1);
        end
        %13/3/2014 modified for B3D multichannel
        if  isempty(mat1)==1,
            d = whos;
            str = {d.name};
            [s,v] = listdlg('PromptString','Select a variable for analysis:',...
                'SelectionMode','single',...
                'ListString',str);
            % la scelta va in s (cell array)
            mat1=eval(char(str(s)));
        end
        
        
        
        % -----------------------------------------------------------------------------------------------   
        
    case 'loadAW',
        [NOME,PERC111]=uigetfile('*.mat','Import File with TW Average to be analyzed');
        if NOME==0,
            warndlg('No File selected !','!! WARNING !!');
        else
            ST0=['PERC','FILE'];
            STR=strrep(ST0,'PERC',PERC111);
            STR1=strrep(STR,'FILE',NOME);
            % Leggo file
            load(STR1);
            q=min(size(Temp));
            template=Temp(q,:);
            % cerco dove finisce
            trovazero=find(template(1,:)<=(template(1,1)-0.5));
            if isempty(trovazero),
                T=length(template);
            else
                T=trovazero(1,1);								%T=campione finale del template
            end;
            T=T-1;
            % modifica introdotta per LAA
            if LAA==1,
                while abs(template(1,T)-template(1,1))>5,
                    T=T-1;
                end;
            end;
            template=template(1,1:T);
            figure(1),plot(template,'r'),ylabel('ml'),xlabel('samples (Fc=400 Hz)'),title('Final Average')
        end
        caricoAW=1;
        % -------------------------------------------------------------------------------------------------   
        
        
    case 'ave',
        if ~isempty(mat1),
            % definisco a priori alcuni parametri, non modificabili dall'utente
            p=0;		% vincolo sulla WF (0 1)
            cif=0;	% 1 per AQ; cifra di merito da utilizzare (0=amp+dx 10=amp 11=dx 1=amp+dx+d2x 2=amp+d2x 3=d2x  4=dx+d2x)
            MMM=0;	% effettua la media sui valori normalizzati e detrendati (ho rimosso il codice che non c'entrava)
            % _______________________________________
            ma2=mat1;
            dimm=size(ma2);
            matorig=mat1;
            clear mat1
            OK=0;
%             keyboard
            while OK==0,
                % inserisci estremi spezzone da analizzare      
                prompt={'Enter the starting cycle number to be analyzed:','Enter the number of cycles to average:'};
                def={'1','32'};
                dlgTitle='Input for Cycle interval to Analyze';
                lineNo=1;
                answer=inputdlg(prompt,dlgTitle,lineNo,def);
                % converto da cell array a char
                answer1=char(answer);
                % controllo sui valori inseriti
                n1=str2num(answer1(1,:));
                N=str2num(answer1(2,:));      
                if (n1>dimm(1,1))|(n1+N-1>dimm(1,1)), %Modificado Alba 27/10/2014
                    errordlg('Wrong selection: those cycles do not exist!','!! WARNING !!');
                else
                    % estraggo i battiti prescelti per calcolare la wf
                    mat1=ma2(n1:(n1+N-1),:);
                    % sposto per poi salvare i riferimenti ai battiti sui quali è stata calcolata la media
                    nAW=n1;
                    NAW=N;
                    % n è il numero di righe di mat1 ed è potenza di 2
                    dimmat=size(mat1);
                    n=dimmat(1,1);
                    M=dimmat(1,2);
                end;
                for i=1:n,
                    fine=find(mat1(i,:)==0);
                    if isempty(fine),
                        fine=M;
                    end
                    plot(mat1(i,1:(fine-1))),clear fine,hold on
                end
                hold off
                ButtonName=questdlg('Do you confirm your selection ?','Confirm selection','Yes');
                if strcmp(ButtonName,'Yes')==1,
                    OK=1;
                end
            end
            % ---------------------
            % inserisci ampiezza finestra di calcolo      
            prompt={'Enter the width of computing window in samples:'};
            def={'10'};
            dlgTitle='Input for window width';
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,def);
            % converto da cell array a char
            answer1=char(answer);
            % controllo sui valori inseriti
            amp=str2num(answer1(1,:));
            if amp>M,
                errordlg('Wrong selection: window width exceed signal length!','!! WARNING !!');
            end
            % ----------------------------------------------------------------------------------
            % PROCEDURA DI CALCOLO DELLA MEDIA NON LINEARE 
            % prendo 2 battiti alla volta e calcolo la wf e la forma d'onda
            % risultante, che viene messa in mat2 
            tot=n-1;													 % totale iterazioni da effettuare	
            itot=1;															 % conta le iterazioni effettuate
            h = waitbar(0,'Please wait...');
            pippo=1;
            while n>1, 														 %n=#righe da analizzare
                waitbar(itot/tot,h)
                % creo la matrice buffer mat2
                mat2=zeros(n/2,M);										%M=#colonne da analizzare
                kappa=1;
                erre=1;														%erre=# di riga esaminata
                for ii=1:n/2
                    s1=0;
                    s2=0;
                    % pongo in s1 la erre-sima riga di mat1 fino a che trova 0 o finisce la riga
                    % verifico di iniziare da un valore diverso da zero
                    jj=1;
                    while mat1(erre,jj)==0,
                        jj=jj+1;												%jj=campione iniziale del ciclo
                    end
                    trovazero=find(mat1(erre,jj:M)==0);
                    if isempty(trovazero),
                        jjj=M;
                    else
                        jjj=trovazero(1,1);								%jjj=campione finale del ciclo
                    end;
                    %17022014: aggiunge uno zero alla fine, lo tolgo con -1
                    s1=mat1(erre,jj:jjj-1);
                    clear jj,clear jjj,clear trovazero
                    
                    % pongo in s2 la (erre+1)-sima riga di mat1 fino a che trova 0 o finisce la riga
                    % verifico di iniziare da un valore diverso da zero
                    jj=1;
                    while mat1(erre+1,jj)==0,
                        jj=jj+1;
                    end
                    trovazero=find(mat1(erre+1,jj:M)==0);
                    if isempty(trovazero),
                        jjj=M;
                    else
                        jjj=trovazero(1,1);								%jjj=campione finale del ciclo
                    end;
                    %17022014: aggiunge uno zero alla fine, lo tolgo con -1
                    s2=mat1(erre+1,jj:jjj-1);
                    
                    % Calcolo warping function con allineamento dei segnali
                    % per 2 battiti x,y
                    % tengo conto della derivata seconda nel calcolo della f.di dissimilarità
                    % versione 13/7/98
                    %---------------------------------------------------------------------------
                    % P0NPROVA
                    % preparazione dei vettori tramite sottrazione della retta
                    % congiungente i punti estremi (per hp. dipendente dal resp)
                    % e normalizzazione tra 0 e 1 delle ampiezze
                    counterr=0;
                    x=s1;
                    y=s2;
                    J=max(size(y));                           
                    I=max(size(x));
                    % creo la retta che congiunge i punti estremi (primo,ultimo)
                    coangx=(x(I)-x(1))/(I-1);
                    quotax=(x(1)*I-x(I))/(I-1);
                    i=linspace(1,I,I);
                    lx=coangx.*i+quotax;
                    coangy=(y(J)-y(1))/(J-1);	
                    quotay=(y(1)*J-y(J))/(J-1);
                    j=linspace(1,J,J);
                    ly=coangy.*j+quotay;
                    clear j,clear i
                    % sottraggo la retta così che i due punti estremi abbiano la stessa ampiezza
                    % 8/9/98 salvo i vettori originari prima del detrend
                    xx=x;														%xx= 1° curva prima del detrend
                    yy=y;     												%yy= 2° curva prima del detrend
                    %17022014: errore in secondo passaggio su curva gia'
                    % detrendata. Retta va a diminuire picchi. Posto
                    % vincolo in funzione coeff ang. nuova retta detrend
                    if abs(coangx)>.0001,
                        x=(x-lx)+min(x);
                    end
                    if abs(coangy)>.0001,
                        y=(y-ly)+min(y);
                    end
                    % 14/1/99
                    % calcolo der.prima sui segnali detrendati e la normalizzo
                    % 5/10/98
                    % derivata prima
                    % creo il vettore dx(i)=x(i)-x(i-1)
                    dx=zeros(1,I);
                    diffx=diff(x);
                    dx(1,2:I)=diffx;
                    % creo il vettore dy(j)=y(j)-y(j-1)
                    dy=zeros(1,J);
                    diffy=diff(y);
                    dy(1,2:J)=diffy;
                    % normalizzo le ampiezze tra 0 e 1
                    % 5/10/98
                    % SE NORMALIZZO QUI TRA 0 E 1 PERDO IL SEGNO!!
                    dx=(dx-min(dx(2:I)))/(max(dx)-min(dx(2:I)));
                    dy=(dy-min(dy(2:J)))/(max(dy)-min(dy(2:J)));
                    
                    % se devo usarla, creo la derivata seconda
                    % 2/10/98
                    if (cif==1)|(cif==2)|(cif==3)|(cif==4),
                        % creo il vettore dx2(i)=x(i)-2*x(i-1)+x(i-2)
                        dx2=zeros(1,I);
                        diffx1=diffx(1,2:length(diffx));
                        dx2(1,3:I)=diffx1-diffx(1,1:length(diffx)-1);                 
                        % creo il vettore dy2(j)=y(j)-2*y(j-1)+y(j-2)
                        dy2=zeros(1,J);
                        diffy1=diffy(1,2:length(diffy));
                        dy2(1,3:J)=diffy1-diffy(1,1:length(diffy)-1);   
                        % normalizzo le ampiezze tra 0 e 1
                        % 5/10/98
                        % SE NORMALIZZO QUI TRA 0 E 1 PERDO IL SEGNO!!
                        dx2=(dx2-min(dx2(3:I)))/(max(dx2)-min(dx2(3:I)));
                        dy2=(dy2-min(dy2(3:J)))/(max(dy2)-min(dy2(3:J)));
                        clear diffx1,clear diffy1
                    end;
                    clear diffx,clear diffy         
                    %17022014: calcolo ofs e delta su segnali originari e
                    %non detrendati (se no sottostima)
                    % definisco i parametri per poter riscalare il template sui segnali detrendati
%                     ofs=mean([min(x),min(y)]);
%                     delta=mean([(max(x)-min(x)),(max(y)-min(y))]);
                    ofs=mean([min(s1),min(s2)]);
                    delta=mean([(max(s1)-min(s1)),(max(s2)-min(s2))]);
                    
                    % normalizzo le ampiezze così che siano tra 0 e 1 
                    x=(x-min(x))/(max(x)-min(x));    
                    y=(y-min(y))/(max(y)-min(y));   
                    clear lx,clear ly
                    % -------------------------------------------------------------------------
                    %P1NPROVA 
                    %P1MOD+DISS
                    % Calcolo la matrice d delle dissimilarità locali dei vettori x,y
                    % per gli elementi all'interno della finestra di ampiezza r
                    % cif=modalità di calcolo della cifra di merito:
                    % 0=amp+dx; 10=amp; 11=dx; 1=amp+dx+d2x; 2=amp+d2x
                    % 3=d2x;     4=dx+d2x;
                    % amp=ampiezza finestra, per segnali di lunghezza uguale, in cui il costo è diverso da 100
                    % va inserito dall'utente.
                    % normalizzo le derivate tra 0 e 1 dopo il calcolo della matrice delle distanze relative
                    % versione 6/10/98
                    
                    if (I>amp)&(J>amp),
                        r=amp+abs(I-J);
                        % 17/2/99
                        if (I<=r),
                            r=I;
                        end
                        if (J<=r),
                            r=J;
                        end;
                    else
                        r=min(I,J);
                    end;
                    % 6/10/98 evito di calcolare d se la cifra di merito non lo richiede
                    if (cif==0)|(cif==10)|(cif==1)|(cif==2),
                        d0=zeros(I,J);
                    end
                    if (cif==0)|(cif==11)|(cif==1)|(cif==4),
                        d1=zeros(I,J);
                    end
                    if (cif>=1)&(cif<=4),
                        d2=zeros(I,J);
                    end
                    % qui inizia il calcolo delle matrici
                    for i=1:I              
                        for j=1:J                     
                            if (i<j+r)&(j<i+r)
                                if (i==1)|(j==1)
                                    if (cif==0)|(cif==10)|(cif==1)|(cif==2),
                                        d0(i,j)=abs(x(i)-y(j));
                                    end;
                                end;
                                if (i>1)&(j>1),
                                    if (cif==0)|(cif==10)|(cif==1)|(cif==2),
                                        d0(i,j)=abs(x(i)-y(j));
                                    end;
                                    if (cif==0)|(cif==11)|(cif==1)|(cif==4),
                                        d1(i,j)=abs(dx(i)-dy(j));
                                    end;
                                    if (i>2)&(j>2)&(cif>=1)&(cif<=4),
                                        d2(i,j)=abs(dx2(i)-dy2(j));
                                    end;
                                end;
                            end;
                        end;
                    end 
                    % costruzione matrice distanza totale come somma matrici ampiezza e derivate
                    D=zeros(I,J);
                    if (cif==0)|(cif==10)|(cif==1)|(cif==2),
                        D=D+d0;
                        clear d0;
                    end;
                    if (cif==0)|(cif==1)|(cif==4)|(cif==11),
                        D=D+d1;
                        clear d1;
                    end;
                    if (cif==1)|(cif==2)|(cif==3)|(cif==4),
                        D=D+d2;
                        clear d2
                    end
                    d=D;
                    clear D
                    %--------------------------------------------------------------------------
                    % P2MOD
                    % inizializzazione della matrice g con vincolo sulla finestra
                    i=1;
                    j=1;
                    g=ones(I,J)*1000;
                    %for i=1:I
                    %   for j=1:J
                    %      if (i>=j+r)|(j>=i+r)
                    %         g(i,j)=1000;
                    %      end; 
                    %  end;
                    %end;
                    %---------------------------------------------------------------------------
                    % P3F
                    % implementazione del metodo di Picton per la creazione
                    % della matrice g con inizializzazione alla Sakoe
                    % "Electroencephal. and Clin. Neuroph., 1988,71:212-225
                    % 23/12/98
                    % rivisto quella per p=1 (qualche problema legato all'inizializzazione)
                    
                    % inizializzazione
                    i=1;
                    j=1;
                    g(1,1)=2*d(1,1);
                    % ----------------------------------------------------
                    % p=0
                    if p==0
                        % cornice orizzontale
                        for j=2:r,
                            g(1,j)=g(1,j-1)+d(1,j);
                        end
                        % cornice verticale
                        for i=2:r,             
                            g(i,1)=g(i-1,1)+d(i,1);
                        end
                        % interno della cornice
                        i=1;j=2;
                        while j<=J,
                            i=i+1;
                            if i>=j+r
                                j=j+1;
                                if j>r
                                    i=j-r;
                                else i=1;
                                end;
                            elseif i<=I
                                g(i,j)=min([g(i-1,j)+d(i,j),
                                    g(i-1,j-1)+2*d(i,j),
                                    g(i,j-1)+d(i,j)]);
                            end; 
                        end;
                        % ------------------------------------------------------
                        % p=1
                    elseif p==1
                        % cornice orizzontale
                        for j=2:r,
                            g(1,j)=g(1,j-1)+d(1,j);
                        end
                        % cornice verticale
                        for i=2:r,
                            g(i,1)=g(i-1,1)+d(i,1);
                        end
                        % 2° cornice verticale
                        j=2;
                        for i=2:r+1,
                            g(i,j)=min([g(i-1,j)+d(i,j),
                                g(i-1,j-1)+2*d(i,j),
                                g(i,j-1)+d(i,j)]);
                        end;
                        % 2° cornice orizzontale
                        i=2;
                        if r>=2,
                            for j=3:r+1,
                                g(i,j)=min([g(i-1,j)+d(i,j),
                                    g(i-1,j-1)+2*d(i,j),
                                    g(i,j-1)+d(i,j)]);
                            end
                        end;
                        % interno della cornice
                        i=2;j=3;
                        while j<=J,
                            i=i+1;
                            if i>=j+r
                                j=j+1;
                                if j>r+2
                                    i=j-r;
                                else i=2;
                                end;
                            elseif i<=I
                                g(i,j)=min([g(i-1,j-2)+2*d(i,j-1)+d(i,j),
                                    g(i-1,j-1)+2*d(i,j),
                                    g(i-2,j-1)+2*d(i-1,j)+d(i,j)]);
                            end; 
                        end;
                        D=g(I,J)/(I+J);
                        %----------------------------------------------------------   
                    end;
                    clear d;
                    %---------------------------------------------------------------------------
                    % PICT4
                    % implementazione del metodo di Picton per la ricerca
                    % della funzione di warping
                    % "Electroencephal. and Clin. Neuroph., 1988,71:212-225
                    
                    % inizializzazione puntatori
                    % la ricerca procede a ritroso da (I,J) a (1,1)
                    i=I;
                    j=J;
                    k=1;
                    wfx(k)=i;
                    wfy(k)=j;
                    k=k+1;
                    while (i>1)&(j>1),
                        gmin=min([g(i-1,j-1),g(i-1,j),g(i,j-1)]);
                        if gmin==g(i-1,j-1),
                            wfx(k)=i-1;
                            wfy(k)=j-1;
                            % 17/2/99
                            if (i-1)<=(j-r)|(j-1)<=(i-r),
                                %warndlg('Error ! Window too narrow!')
                                counterr=counterr+1;
                            end;
                            i=i-1;
                            j=j-1;
                        elseif gmin==g(i-1,j),
                            wfx(k)=i-1;
                            if (i-1)<=(j-r),
                                %warndlg('Error ! Window too narrow!')
                                counterr=counterr+1;
                            end;
                            wfy(k)=j;
                            i=i-1;
                        elseif gmin==g(i,j-1),
                            wfx(k)=i;
                            wfy(k)=j-1;
                            if (j-1)<=(i-r),
                                %warndlg('Error ! Window too narrow!')
                                counterr=counterr+1;
                            end;
                            j=j-1;
                        end;
                        k=k+1;
                    end;
                    % sono sulla cornice orizzontale ma non nell'origine
                    if (i==1)&(j>1),
                        while j>1,
                            wfx(k)=i;
                            wfy(k)=j-1;
                            j=j-1;
                            k=k+1;
                        end;
                        % sono sulla cornice verticale ma non nell'origine
                    elseif (i>1)&(j==1),
                        while i>1,
                            wfx(k)=i-1;
                            wfy(k)=j;
                            i=i-1;
                            k=k+1;
                        end;
                    end;
                    clear g;
                    
                    %________________________________________________
                    % SE VOGLIO VISUALIZZARE DEVO TOGLIERE I COMMENTI
%                     
%                     keyboard
%                     figure
%                     plot(wfx,wfy)
% %                  %  traccio i punti omologhi
%                     K=max(size(wfx));
%                     points=K;
%                     figure(2)
%                     plot(x),hold on
%                     plot(y,'r')
%     %                plot(1:max(size(s1)),s1(1:max(size(s1))),1:max(size(s2)),s2(1:max(size(s2))),'r'),hold on
%                     for i=1:points,
%                     	T1=[wfx(i) wfy(i)];
%                     	T2=[x(wfx(i)) y(wfy(i))];
%                        % T2=[s1(wfx(i)) s2(wfy(i))];
%                     	line(T1,T2)
%                     	hold on
%                     end
%                     pause
%                     hold off
%                     points=K;
%                     figure(3)
%                     dif1=[0 diff(s1)];
%                     dif2=[0 diff(s2)];
%                     plot(1:max(size(dif1)),dif1(1:max(size(dif1))),1:max(size(dif2)),dif2(1:max(size(dif2)))),hold on
%                     for i=1:points,
%                     	T1=[wfx(i) wfy(i)];
%                     	T2=[dif1(wfx(i)) dif2(wfy(i))];
%                     	line(T1,T2)
%                     	hold on
%                     end
%                     pause
%                     hold off
                    %-------------------------------------------------------------------------
                    
                    % PICT5
                    % implementazione del metodo per la costruzione
                    % del template tra i due vettori a partire dalla wf calcolata
                    % e successiva riscalatura delle ampiezze ai valori originari
                    
                    % il vettore v conterrà gli elementi dei due vettori di partenza 
                    % combinati secondo la wf:
                    % - se ho uno spostamento orizzontale o verticale prendo x(i),y(j)
                    % - se ho uno spostamento diagonale prendo x(i),y(j),x(i),y(j)
                    
                    K=max(size(wfx));
                    v(1)=x(wfx(K));
                    v(2)=y(wfy(K));
                    % 4/9/98 voglio che il primo elemento sia la media dei due primi elementi
                    v(3)=x(wfx(K));
                    v(4)=y(wfy(K));
                    s=5;
                    k=K;
                    while k>1,
                        % spostamento diagonale 
                        if abs(wfx(k)-wfx(k-1))==abs(wfy(k)-wfy(k-1))
                            for q=1:2,
                                v(s)=x(wfx(k-1));
                                s=s+1;
                                v(s)=y(wfy(k-1));
                                s=s+1;
                            end;
                            % spostamento orizzontale
                        elseif abs(wfx(k)-wfx(k-1))>abs(wfy(k)-wfy(k-1))
                            v(s)=x(wfx(k-1));
                            s=s+1;
                            v(s)=y(wfy(k-1));
                            s=s+1;
                            % spostamento verticale
                        elseif abs(wfx(k)-wfx(k-1))<abs(wfy(k)-wfy(k-1))
                            v(s)=x(wfx(k-1));
                            s=s+1;
                            v(s)=y(wfy(k-1));
                            s=s+1; 
                        end;
                        k=k-1;
                    end;
                    
                    % calcolo media per la prima volta
                    j=1;
                    for i=1:max(size(v))/2
                        A(i)=mean([v(j),v(j+1)]);
                        j=j+2;
                    end;
                    % calcolo media per la seconda volta
                    % 4/9/98 vedi sopra
                    dimA=max(size(A));
                    % se il numero di elementi di A è dispari, medio dal 2°elemento di A
                    if rem(dimA,2)>0,
                        B(1)=A(1);
                        j=2;
                        for i=2:(dimA+1)/2
                            B(i)=mean([A(j),A(j+1)]);
                            j=j+2;
                        end;
                        
                    else
                        j=1;
                        for i=1:dimA/2
                            B(i)=mean([A(j),A(j+1)]);
                            j=j+2;
                        end;
                    end
                    clear A;
%                     % 17022014: visualizzo media risultante prima del
%                     % detrend
%                     figure
%                     plot(x),hold on
%                     plot(y,'b')
%                     plot(B,'r')
%                     hold off
                    % riscalo il template trovato con i valori medi del minimo
                    % e del (max-min) calcolati da x e y dopo il detrend lineare
                    a=B.*delta+ofs;
                    
                    %17022014: se detrendo trovo onda non piú intermedia
                    % per cui tolgo questa parte. Media su segnali
                    % detrendati
%                     % 24/9/98
%                     % riapplico il detrend mediando i coefficenti dei due segnali originari
%                     coang=(coangx+coangy)/2;
%                     quota=(quotax+quotay)/2;
%                     i=linspace(1,max(size(B)),max(size(B)));
%                     lx=coang.*i+quota;
%                     % sommo la retta così da riottenere per i due punti estremi
%                     % l'ampiezza mediata corrispondente
%                     a=(a+lx)-min(a);
                    A=max(size(a));
                    % pongo la forma d'onda generata nella kappa-sima riga di mat2 
                    mat2(kappa,1:A)=a;
                   
                    
                    % passo ai successivi 2 battiti
                    kappa=kappa+1;
                    erre=erre+2;
                    itot=itot+1;
                    waitbar(itot/tot,h)
                    clear a,clear v,clear wfx,clear wfy,clear dx,clear dy,clear dx2,
                    clear dy2,clear xdetr,clear ydetr,clear x,clear y,clear s1,clear s2,
                    clear B, clear xx,clear yy,clear coangx,clear coangy
                    clear quotax,clear quotay,clear quota,clear coang,clear lx
                end; 
                mat1=mat2;
                
                % copio in Temp la prima riga di mat2
                % per estrarre i template ottenuti con 2^pippo battiti
                Temp(pippo,1:M)=mat2(1,1:M);
                pippo=pippo+1;
                clear mat2;
                % aggiorno il numero di righe di mat1 (e mat2)
                n=size(mat1);
                n=n(1);   
            end;
            drawnow;close(h)
            %Visualizzazione Template finale estratto
            %dimtemp=size(Temp);
            q=min(size(Temp));
            template=Temp(q,:);
            %17022014: visualizzo solo non zeri
            fint=find(template);
            sfint=length(fint);
            figure(1),plot(template(1:sfint),'r'),ylabel('ml'),xlabel('samples'),title('Final Average') 
            counterr
            clear counterr
            % abilito calcolo dei parametri
            caricoAW=1; carico=1; %060618 Federica, per calolare TWF
            % salva matrice Temp
            ButtonName=questdlg('Save results on .mat file?','SAVE in Matlab format','Yes');
            if strcmp(ButtonName,'Yes')==1,
                
                save('DTWresults.mat');
            end;
        else
            errordlg('No file imported in work area!','!!! WARNING !!!');
        end;
        
        % --------------------------------------------------------------------------------------------------   
    case 'extract'
        if ~isempty(Temp),
            q=min(size(Temp));
            qq=max(size(Temp));
            template=Temp(q,:);
            % cerco dove finisce
            trovazero=find(template(1,:)<=(template(1,1)-0.5));
            if isempty(trovazero),
                T=length(template);
            else
                T=trovazero(1,1);		
                % 6/7/01 introdotto modifica per evitare zeri all'inizio
                tz=1;
                while T<100,
                    tz=tz+1;
                    T=trovazero(1,tz);
                end
                
                %T=campione finale del template
            end;
            T=T-1;
            % modifica introdotta per LAA
            if LAA==1,
                while abs(template(1,T)-template(1,1))>5,
                    T=T-1;
                end;
            end;
            template=template(1,1:T);
            figure(1),plot(template,'r'),ylabel('ml'),xlabel('samples (Fc=300 Hz)'),title('Final Average')
            
            % Frequenza di campionamento del segnale AQ (in Hz)
            % cerco dove finisce
            % leggo il filtro derivativo
            load c:\fd10_300;
            temp=[template(1:T) template(1:T) template(1:T) template(1:T)];
            dtemp=filter(coeff,1,temp);
            
            % cerco il max corrispondente all'EDV
            edv=max(template);
            [t1,tedv]=find(template==edv);
            tedv=tedv(1,1);
            % estraggo la derivata corretta buttando via i primi 75 campioni e prendendo il template centrale
            dtemplate=dtemp(76+T:76+2*T); 
            % diamo un occhio alla derivata seconda
            ddtemp=[dtemplate dtemplate dtemplate];
            dtemp2=filter(coeff,1,ddtemp);
            dtemplate2=dtemp2(76+T:76+2*T); 
            
            % ricerca dei punti fiduciari sul template finale 
            edv;	% ml
            esv=template(1);	% ml
            sv=edv-template(T);	% ml
            fe=(sv/edv)*100; % a
            pfr=max(dtemplate(1:round(.25*T)));	% ml/sec
            ivrt=min(dtemplate(1:round(.25*T)));	% ml/sec
            %knee=min(dtemplate(round(.25*T):round(.4*T)));	
            pafr=max(dtemplate(round(.4*T):T));	% ml/sec
            per=min(dtemplate(round(.4*T):round(.8*T)));	% ml/sec
            indvol=[edv esv sv fe]
            
            % ricerca degli indici temporali
            % la scala dei tempi parte dall'ESV ed è trasformata in msec
            tedv=tedv/Fc;
            tknee=1;
            [b,c]=find(dtemplate==pfr);
            tpfr=round(c(1,1))/Fc;	% sec
            clear b,clear c
            % ____________________________
            % IVRT
            % inizio IVRT
            tivrts=1/Fc;		%sec
            % IVRT fine
            [b,c]=find(dtemplate==ivrt);
            tivrte=round(c(1,1))/Fc;	% sec
            clear b,clear c
            
            %tknee=find(dtemplate(round(.25*T):round(.4*T))==knee);
            %tknee=tknee(1,1)/Fc;	% sec
            [b,c]=find(dtemplate==pafr);
            tpafr=round(c(1,1))/Fc;	% sec
            clear b,clear c
            [b,c]=find(dtemplate==per);
            tper=round(c(1,1))/Fc;	% sec
            clear b,clear c
            % calcolo terzo TA
            diastole=template(1:(tedv*Fc))-esv;
            centoT=tedv*Fc;
            centoA=edv-esv;
            terzoT=round(centoT/3);
            terzoA=centoA/3;
            time13=100*(diastole(terzoT)/centoA)
            difA=abs(diastole-terzoA);mi=min(difA);
            cerca=find((diastole(:)==mi+terzoA)|(diastole(:)==terzoA-mi));
            amp13=100*cerca/(tedv*Fc)
            % visualizzo su figura 2
            figure(2)
            set(gcf,'PaperType','a4letter');
            pupu=100*(diastole(1:(tedv*Fc)))/sv;
            plot(pupu);axis([0 tedv*Fc 0 100]);
            xlabel('% Diastolic Filling period'),ylabel('% Diastolic Filling')
            hold on
            x1=[terzoT terzoT];
            y1=[0 time13];line(x1,y1);
            x2=[0 cerca];
            y2=[100*terzoA/sv 100*terzoA/sv];
            line(x2,y2)
            hold off
            
            % cerco il minimo della derivata seconda tra pfr e knee
            % con Wk definisco l'ampiezza in campioni della finestra, a partire da pfr,
            % nella quale procedo al calcolo del minimo della derivata seconda
            Wk=100;
            kneeD2=min(dtemplate2(round(tpfr*Fc):round(tpfr*Fc)+Wk));
            [b,c]=find(dtemplate2==kneeD2);
            tkneeD2=round(c(1,1))/Fc	% sec
            clear b, clear c
            
            % cerco il massimo della derivata seconda tra knee2 e pafr
            knee2D2=max(dtemplate2(round(tkneeD2*Fc):round(tpafr*Fc)));
            if isempty(knee2D2),
                knee2D2=dtemplate2(tpafr*Fc);
            end;
            [b,c]=find(dtemplate2==knee2D2);
            tknee2D2=round(c(1,1))/Fc;	
            
            % visualizzo i punti di knee individuati
            f=[template(round(tkneeD2*Fc)) template(round(tknee2D2*Fc))];
            e=[round(tkneeD2*Fc) round(tknee2D2*Fc)];
            g=[dtemplate2(round(tkneeD2*Fc)) dtemplate2(round(tknee2D2*Fc))];
            d=[dtemplate(round(tivrts*Fc)) dtemplate(round(tivrte*Fc)) dtemplate(round(tpfr*Fc)) dtemplate(round(tpafr*Fc)) dtemplate(round(tper*Fc))];
            c=[round(tivrts*Fc) round(tivrte*Fc) round(tpfr*Fc) round(tpafr*Fc) round(tper*Fc)];
            b=[template(round(tivrts*Fc)) template(round(tivrte*Fc)) template(round(tpfr*Fc)) template(round(tpafr*Fc)) template(round(tper*Fc))];
            figure(3);
            set(gcf,'PaperType','a4letter');
            subplot(311),plot(1:T,template(1:T)),title('template'),ylabel('ml'),grid on
            hold on
            subplot(311),plot(e,f,'mx'),plot(c,b,'go');
            hold off
            subplot(312),plot(1:T,dtemplate(1:T)),title('first derivative'),ylabel('ml/sec'),xlabel('# sample'),grid on;
            hold on
            subplot(312),plot(c,d,'go');
            hold off
            subplot(313),plot(1:T,dtemplate2(1:T)),title('second derivative'),ylabel('ml/sec^2'),xlabel('# sample'),grid on;
            hold on
            subplot(313),plot(e,g,'mx');
            hold off
            
            % Richiesta per correzione
            OK=0;
            while OK==0,
                ButtonName=questdlg('Would you correct any point?','Correction','No');
                if strcmp(ButtonName,'Yes')==1,
                    % Correzione punti individuati      
                    prompt={'Enter the kind of correction: knee(1), knee2(2), both(3), pafr(4), pafr and knee1(5), pafr and knee2(6), pfr(7), per(8), IVRTstart(9), IVRTend(10)'};
                    def={''};
                    dlgTitle='Input for correction';
                    lineNo=1;
                    answer=inputdlg(prompt,dlgTitle,lineNo,def);
                    % converto da cell array a char
                    answer1=char(answer);
                    % controllo sui valori inseriti
                    MOD=str2num(answer1(1,:));
                    if (MOD>11)|(MOD<1),
                        errordlg('Wrong selection!','!! WARNING !!');
                    end
                    subplot(111);
                    % visualizzo sovrapposti il segnale di volume e le sue derivate
                    % Figura4 
                    if (MOD==1)|(MOD==2)|(MOD==3)|(MOD==4)|(MOD==5)|(MOD==6),
                        figure(4)
                        plot(1:T,template(1:T)),grid on;
                        hold on
                        plot(c,b,'mx');
                        hold on
                        plot(e,f,'go');
                        hold off
                    end;
                    if (MOD==1)|(MOD==3)|(MOD==5),
                        % inserisci il valore di knee tramite mouse
                        figure(4)
                        helpdlg('Press Enter and choose the knee point','Knee Point Selection');pause
                        [tkneeD2,kneeD2]=ginput(1);
                        tfpfr=round(tkneeD2);
                        tkneeD2=round(tkneeD2)/Fc;	%sec
                        %kneeD2=template(round(tkneeD2*Fc));
                        f=[template(round(tkneeD2*Fc)) template(round(tknee2D2*Fc))];
                        e=[round(tkneeD2*Fc) round(tknee2D2*Fc)];
                        g=[dtemplate2(round(tkneeD2*Fc)) dtemplate2(round(tknee2D2*Fc))];
                        figure(3);
                        subplot(311),plot(1:T,template(1:T)),title('template'),ylabel('ml'),grid on
                        hold on
                        subplot(311),plot(e,f,'mx'),plot(c,b,'go');
                        hold off
                        subplot(313),plot(1:T,dtemplate2(1:T)),title('second derivative'),ylabel('ml/sec^2'),xlabel('# sample'),grid on;
                        hold on
                        subplot(313),plot(e,g,'mx');
                        hold off     
                    end
                    
                    if (MOD==2)|(MOD==3)|(MOD==6),
                        % inserisci il valore di knee2 tramite mouse
                        figure(4)
                        helpdlg('Press Enter and choose the knee2 point','Knee2 Point Selection');pause
                        [tknee2D2,knee2D2]=ginput(1);
                        tipafr=round(tknee2D2);
                        tknee2D2=round(tknee2D2)/Fc;	%sec
                        %knee2D2=template(round(tknee2D2*Fc));
                        f=[template(round(tkneeD2*Fc)) template(round(tknee2D2*Fc))];
                        e=[round(tkneeD2*Fc) round(tknee2D2*Fc)];
                        g=[dtemplate2(round(tkneeD2*Fc)) dtemplate2(round(tknee2D2*Fc))];
                        figure(3);
                        subplot(311),plot(1:T,template(1:T)),title('template'),ylabel('ml'),grid on
                        hold on
                        subplot(311),plot(e,f,'mx'),plot(c,b,'go');
                        hold off
                        subplot(313),plot(1:T,dtemplate2(1:T)),title('second derivative'),ylabel('ml/sec^2'),xlabel('# sample'),grid on;
                        hold on
                        subplot(313),plot(e,g,'mx');
                        hold off               
                    end
                    
                    if (MOD==4)|(MOD==5)|(MOD==6),
                        % inserisci il valore di pafr tramite mouse
                        figure(4)
                        plot(1:T,dtemplate(1:T))
                        helpdlg('Press Enter and choose the pafr point','Peak atrial Filling Rate Point Selection');pause
                        [tpafr,pafr]=ginput(1);
                        % in msec
                        tpafr=round(tpafr)/Fc;	%sec
                        pafr=dtemplate(round(tpafr*Fc));
                        d=[dtemplate(round(tpfr*Fc)) dtemplate(round(tpafr*Fc)) dtemplate(round(tper*Fc))];
                        c=[round(tpfr*Fc) round(tpafr*Fc) round(tper*Fc)];
                        b=[template(round(tpfr*Fc)) template(round(tpafr*Fc)) template(round(tper*Fc))];
                        figure(3);
                        subplot(311),plot(1:T,template(1:T)),title('template'),ylabel('ml'),grid on
                        hold on
                        subplot(311),plot(e,f,'mx'),plot(c,b,'go');
                        hold off
                        subplot(312),plot(1:T,dtemplate(1:T)),title('first derivative'),ylabel('ml/sec'),xlabel('# sample'),grid on;
                        hold on
                        subplot(312),plot(c,d,'go');
                        hold off                          
                    end
                    
                    if MOD==7,
                        % inserisci il valore di pfr tramite mouse
                        figure(4)
                        plot(1:T,dtemplate(1:T))
                        helpdlg('Press Enter and choose the pfr point','Peak Filling Rate Point Selection');pause
                        [tpfr,pfr]=ginput(1);
                        % in msec
                        tpfr=round(tpfr)/Fc;	%sec
                        pfr=dtemplate(round(tpfr*Fc));
                        d=[dtemplate(round(tpfr*Fc)) dtemplate(round(tpafr*Fc)) dtemplate(round(tper*Fc))];
                        c=[round(tpfr*Fc) round(tpafr*Fc) round(tper*Fc)];
                        b=[template(round(tpfr*Fc)) template(round(tpafr*Fc)) template(round(tper*Fc))];
                        figure(3);
                        subplot(311),plot(1:T,template(1:T)),title('template'),ylabel('ml'),grid on
                        hold on
                        subplot(311),plot(e,f,'mx'),plot(c,b,'go');
                        hold off
                        subplot(312),plot(1:T,dtemplate(1:T)),title('first derivative'),ylabel('ml/sec'),xlabel('# sample'),grid on;
                        hold on
                        subplot(312),plot(c,d,'go');
                        hold off                          
                    end
                    if MOD==8,
                        % inserisci il valore di per tramite mouse
                        figure(4)
                        plot(1:T,dtemplate(1:T))
                        helpdlg('Press Enter and choose the per point','Peak Ejection Rate Point Selection');pause
                        [tper,per]=ginput(1);
                        % in msec
                        tper=round(tper)/Fc;	%sec
                        per=dtemplate(round(tper*Fc));
                        d=[dtemplate(round(tpfr*Fc)) dtemplate(round(tpafr*Fc)) dtemplate(round(tper*Fc))];
                        c=[round(tpfr*Fc) round(tpafr*Fc) round(tper*Fc)];
                        b=[template(round(tpfr*Fc)) template(round(tpafr*Fc)) template(round(tper*Fc))];
                        figure(3);
                        subplot(311),plot(1:T,template(1:T)),title('template'),ylabel('ml'),grid on
                        hold on
                        subplot(311),plot(e,f,'mx'),plot(c,b,'go');
                        hold off
                        subplot(312),plot(1:T,dtemplate(1:T)),title('first derivative'),ylabel('ml/sec'),xlabel('# sample'),grid on;
                        hold on
                        subplot(312),plot(c,d,'go');
                        hold off                          
                    end
                    if MOD==9,
                        % inserisci il valore di IVRTstart tramite mouse
                        figure(4)
                        plot(1:T,dtemplate(1:T))
                        helpdlg('Press Enter and choose the IVRTstart point','IVRT start Point Selection');pause
                        [tivrts,ivrts]=ginput(1);
                        % in msec
                        tivrts=round(tivrts)/Fc;	% sec
                        % aggiorno automaticamente la fine
                        tivrte=T/Fc;		% sec
                        d=[dtemplate(round(tivrts*Fc)) dtemplate(round(tivrte*Fc)) dtemplate(round(tpfr*Fc)) dtemplate(round(tpafr*Fc)) dtemplate(round(tper*Fc))];
                        c=[round(tivrts*Fc) round(tivrte*Fc) round(tpfr*Fc) round(tpafr*Fc) round(tper*Fc)];
                        b=[template(round(tivrts*Fc)) template(round(tivrte*Fc)) template(round(tpfr*Fc)) template(round(tpafr*Fc)) template(round(tper*Fc))];
                        figure(3);
                        subplot(311),plot(1:T,template(1:T)),title('template'),ylabel('ml'),grid on
                        hold on
                        subplot(311),plot(e,f,'mx'),plot(c,b,'go');
                        hold off
                        subplot(312),plot(1:T,dtemplate(1:T)),title('first derivative'),ylabel('ml/sec'),xlabel('# sample'),grid on;
                        hold on
                        subplot(312),plot(c,d,'go');
                        hold off                          
                    end
                    if MOD==10,
                        % inserisci il valore di IVRTend tramite mouse
                        figure(4)
                        plot(1:T,dtemplate(1:T))
                        helpdlg('Press Enter and choose the IVRTend point','IVRT end Point Selection');pause
                        [tivrte,ivrte]=ginput(1);
                        % in msec
                        tivrte=round(tivrte)/Fc;	% sec
                        % aggiorno automaticamente l'inizio
                        tivrts=1/Fc;	% sec
                        d=[dtemplate(round(tivrts*Fc)) dtemplate(round(tivrte*Fc)) dtemplate(round(tpfr*Fc)) dtemplate(round(tpafr*Fc)) dtemplate(round(tper*Fc))];
                        c=[round(tivrts*Fc) round(tivrte*Fc) round(tpfr*Fc) round(tpafr*Fc) round(tper*Fc)];
                        b=[template(round(tivrts*Fc)) template(round(tivrte*Fc)) template(round(tpfr*Fc)) template(round(tpafr*Fc)) template(round(tper*Fc))];
                        figure(3);
                        subplot(311),plot(1:T,template(1:T)),title('template'),ylabel('ml'),grid on
                        hold on
                        subplot(311),plot(e,f,'mx'),plot(c,b,'go');
                        hold off
                        subplot(312),plot(1:T,dtemplate(1:T)),title('first derivative'),ylabel('ml/sec'),xlabel('# sample'),grid on;
                        hold on
                        subplot(312),plot(c,d,'go');
                        hold off                          
                    end
                    
                else
                    OK=1;
                end
                
            end;
            d=[dtemplate(round(tivrts*Fc)) dtemplate(round(tivrte*Fc)) dtemplate(round(tpfr*Fc)) dtemplate(round(tpafr*Fc)) dtemplate(round(tper*Fc))];
            c=[round(tivrts*Fc) round(tivrte*Fc) round(tpfr*Fc) round(tpafr*Fc) round(tper*Fc)];
            b=[template(round(tivrts*Fc)) template(round(tivrte*Fc)) template(round(tpfr*Fc)) template(round(tpafr*Fc)) template(round(tper*Fc))];
            set(gcf,'PaperType','a4letter','NumberTitle','Off');
            figure('Name',NOME);
            subplot(211),plot(1:T,template(1:T)),title(NOME),ylabel('ml'),grid on
            hold on
            subplot(211),plot(e,f,'mx'),plot(c,b,'go');
            hold off
            subplot(212),plot(1:T,dtemplate(1:T)),title('first derivative'),ylabel('ml/sec'),xlabel('# sample'),grid on;
            hold on
            subplot(212),plot(c,d,'go');
            hold off 
            % stampa
            ButtonName=questdlg('Would you print this figure?','Correction','Yes');
            if strcmp(ButtonName,'Yes')==1,
                set(gcf,'PaperType','a4letter');
                print
            end
            % 23/05/2000
            % durata IVRT
            Tivrt=tivrte-tivrts;		% sec
            if tivrts==1/Fc,
                tpfr=tpfr-tivrte;		% sec
            end
            
            kneeD2=template(round(tkneeD2*Fc));
            knee2D2=template(round(tknee2D2*Fc));
            indder=[pfr pafr per]
            atrfill=edv-template(tknee2D2*Fc);	% ml
            RF100=(kneeD2-esv)/(edv-esv);	% a
            DS100=(knee2D2-kneeD2)/(edv-esv); 	% a
            AF100=atrfill/(edv-esv);		% a
            indfasi=[kneeD2 knee2D2 RF100 DS100 AF100 atrfill]
            tempider=[tpfr tpafr tper]
            ciclo=T*(1/Fc);		% sec
            
            % 23/05/2000
            if tivrts==1/Fc,
                tdiastole=tedv;		% sec
            else
                tdiastole=tedv+Tivrt;	%sec
            end;
            
            tsistole=ciclo-tdiastole;	% sec
            tdiastasi=tknee2D2-tkneeD2;	% sec
            tatrfil=tdiastole-tdiastasi-tkneeD2;
            
            % 23/05/2000
            tPER=tper-tedv;			
            
            tPAFR=tpafr-(tkneeD2+tdiastasi);
            tempifasi=[tkneeD2 tknee2D2 ciclo tdiastole tsistole tdiastasi tatrfil tPER tPAFR Tivrt]
            % abilito il calcolo delle WF tramite TWDYN0
            caricoPAR=1;
            %helpdlg('Results are showed in the Matlab Command Window','RESULTS'),pause 
            
            % visualizza sul grafico i punti fiduciari individuati
            % b=[edv template(round(tpfr*Fc)) template(round(tpafr*Fc)) template(round(tper*Fc))];
            % a=[round(tedv*Fc) round(tpfr*Fc) round(tpafr*Fc) round(tper*Fc)];
            % d=[dtemplate(tedv*Fc) dtemplate(tpfr*Fc) dtemplate(tpafr*Fc) dtemplate(tper*Fc)];
            % c=[round(tedv*Fc) round(tpfr*Fc) round(tpafr*Fc) round(tper*Fc)];
            % f=[template(round(tkneeD2*Fc)) template(round(tknee2D2*Fc))];
            % e=[round(tkneeD2*Fc) round(tknee2D2*Fc)];
            % g=[dtemplate(round(tkneeD2*Fc)) dtemplate(round(tknee2D2*Fc))];
            % set(gcf,'PaperType','a4letter');
            % figure(3);
            % set(gcf,'PaperType','a4letter');
            % subplot(211),plot(1:T,template(1:T)),title('template'),ylabel('ml'),grid on
            % hold on
            % subplot(211),plot(a,b,'mo');
            % hold on
            % subplot(211),plot(e,f,'mx');
            % subplot(212),plot(1:T,dtemplate(1:T)),title('first derivative'),ylabel('ml/sec'),xlabel('# sample'),grid on;
            % hold on
            % subplot(212),plot(c,d,'mo');
            % hold on
            % subplot(212),plot(e,g,'mx');
            % hold off
            
            % visualizzazione derivata seconda
            % subplot(313),plot(1:T,dtemplate2(1:T)),title('second derivative'),ylabel('ml/sec^2'),xlabel('# sample'),grid on;
            % hold on
            % subplot(313),plot(e,f,'mo');
            % hold off
            
            % richiedo all'utente se vuole salvare i risultati
            ButtonName=questdlg('Save results on .xls file?','SAVE in Excel format','Yes');
            if strcmp(ButtonName,'Yes')==1,
                if newfile==0,
                    warndlg('No File selected !','!! WARNING !!');
                end;
                percorso='perc\nome.estensione';
                percorso=strrep(percorso,'nome',newfile);
                percorso=strrep(percorso,'perc',newpath);
                percorso_xls=strrep(percorso,'estensione','xls');
                fid=fopen(percorso_xls,'a');		;
                
                titolo1='ESV(ml)'; titolo2='EDV(ml)'; titolo3='SV(ml)'; titolo4='EF%'; titolo5='KNEE(ml)'; titolo6='KNEE2(ml)';titolo7='PFR(ml/sec)'; titolo8='PAFR(ml/sec)'; titolo9='PER(ml/sec)';titolo10='RF%'; titolo11='DS%'; titolo12='AF%';titolo13='AF(ml)'; titolo14='TPFR(msec)'; titolo15='TPAFR(msec)';titolo16='TPER(msec)'; titolo17='TKNEE(msec)'; titolo18='TKNEE2(msec)';titolo19='CICLO(msec)'; titolo20='Diastole(msec)'; titolo21='Systole(msec)';titolo22='Diastasis(msec)'; titolo23='1/3Time(%Filling)'; titolo24='F1/3FF(%tDiastole)';titolo25='tAF(msec)';titolo26='tPER(msec)';titolo27='tPAFR(msec)';titolo28='Tivrt';
                fprintf(fid,'%s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t %s\t',titolo1,titolo2,titolo3,titolo4,titolo5,titolo6,titolo7,titolo8,titolo9,titolo10,titolo11,titolo12,titolo13,titolo14,titolo15,titolo16,titolo17,titolo18,titolo19,titolo20,titolo21,titolo22,titolo23,titolo24,titolo25,titolo26,titolo27,titolo28);
                fprintf(fid,'\n');
                fprintf(fid,'%5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f\t %5.2f',edv,esv,sv,fe,kneeD2,knee2D2,pfr,pafr,per,RF100*100,DS100*100,AF100*100,atrfill,tpfr*1000,tpafr*1000,tper*1000,tkneeD2*1000,tknee2D2*1000,ciclo*1000,tdiastole*1000,tsistole*1000,tdiastasi*1000,time13,amp13,tatrfil*1000,tPER*1000,tPAFR*1000,Tivrt*1000);
                fprintf(fid,'\n');
                % salvo i segnali per visualizzazione
                titolo1='Amplitude';titolo2='First derivative';
                fprintf(fid,'%s\t',titolo1,titolo2);
                fprintf(fid,'\n');
                for i=1:T,
                    fprintf(fid,'%5.2f\t %5.2f\t',template(i),dtemplate(i));
                    fprintf(fid,'\n');
                end;
                fclose(fid);
                % salvo i riferimenti ai battiti per i quali è stata calcolata la media
                FILEX=['punLV' newfile];
                save('DTWresults.mat');
            end;
            
        else
            errordlg('No average imported in work area!','!!! WARNING !!!');
        end;
        % ---------------------------------------------------------------------------------
    case 'computeTW'
        if (carico==1)&(caricoAW==1),
            % definisco a priori alcuni parametri, non modificabili dall'utente
            dimm=size(matorig);
            NORMA=1;
            p=0;		% vincolo sulla WF (0 1)
            cif=1;	% cifra di merito da utilizzare (0=amp+dx 10=amp 11=dx 1=amp+dx+d2x 2=amp+d2x 3=d2x  4=dx+d2x)
            MMM=0;	% effettua la media sui valori normalizzati e detrendati (ho rimosso il codice che non c'entrava)
            % inserisci ampiezza finestra di calcolo      
            prompt={'Enter the width of computing window in samples:'};
            def={'10'};
            dlgTitle='Input for window width';
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,def);
            % converto da cell array a char
            answer1=char(answer);
            % controllo sui valori inseriti
            amp=str2num(answer1(1,:));
            if amp>M,
                errordlg('Wrong selection: window width exceed signal length!','!! WARNING !!');
            end
            % inserisci estremi spezzone da analizzare      
            prompt={'Enter the starting cycle number to be analyzed:','Enter the final cycle number:'};
            a=num2str(dimm(1,1));
            def={'1',a};
            dlgTitle='Input for Cycle interval to Analyze';
            lineNo=1;
            answer=inputdlg(prompt,dlgTitle,lineNo,def);
            % converto da cell array a char
            answer1=char(answer);
            % controllo sui valori inseriti
            nWF=str2num(answer1(1,:));
            NWF=str2num(answer1(2,:));      
            if (NWF>dimm(1,2))|(nWF>dimm(1,2)),
                errordlg('Wrong selection: those cycles do not exist!','!! WARNING !!');
            else
                % estraggo i battiti prescelti per calcolare la wf
                mat2=matorig(nWF:(nWF+NWF-1),:);
                % n è il numero di righe di mat2 
                dimmat=size(mat2);
                n2=dimmat(1,1);
                M2=dimmat(1,2);
            end;
            
            % Inizializzo
            WFX=zeros(n2,M2);WFY=zeros(n2,M2);
            s1=template;      
            I=max(size(s1));
            % calcolo WF sui valori con detrend e normalizzazione (NORMA=1)
            if NORMA==1,
                % preparazione del vettore template tramite sottrazione della retta
                % congiungente i punti estremi (per hp. dipendente dal resp)
                % e normalizzazione tra 0 e 1 delle ampiezze                
                % creo la retta che congiunge i punti estremi (primo,ultimo)
                coangx=(s1(I)-s1(1))/(I-1);
                quotax=(s1(1)*I-s1(I))/(I-1);
                i=linspace(1,I,I);
                lx=coangx.*i+quotax;         
                % sottraggo la retta così che i due punti estremi abbiano
                % la stessa ampiezza
                s1=(s1-lx)+max(s1);
            end
            clear lx,clear quotax,clear coangx,clear i
            
            % 5/10/98
            % derivata prima
            % se devo usarla creo il vettore dx(i)=x(i)-x(i-1)
            x=s1;        
            if (cif==0)|(cif==11)|(cif==1)|(cif==4),
                dx=zeros(1,I);
                diffx=diff(x);
                dx(1,2:I)=diffx;
                % normalizzo la derivata prima tra 0 e 1
                dx=(dx-min(dx(2:I)))/(max(dx)-min(dx(2:I)));
            end;
            % se devo usarla, creo la derivata seconda
            if (cif==1)|(cif==2)|(cif==3)|(cif==4),
                % creo il vettore dx2(i)=x(i)-2*x(i-1)+x(i-2)
                dx2=zeros(1,I);
                diffx1=diffx(1,2:length(diffx));
                dx2(1,3:I)=diffx1-diffx(1,1:length(diffx)-1);   
                % normalizzo la derivata seconda tra 0 e 1
                % 5/10/98
                % SE NORMALIZZO QUI TRA 0 E 1 PERDO IL SEGNO!!
                dx2=(dx2-min(dx2(3:I)))/(max(dx2)-min(dx2(3:I)));
                clear diffx1,clear diffx
            end;
            % normalizzo le ampiezze così che siano tra 0 e 1 
            x=(x-min(x))/(max(x)-min(x));    
            clear lx,clear coangx,clear quotax
            % inizializzo vettore distanze
            DI=0;
            % pongo in s2 un battito alla volta
            tot=n2;													 % totale iterazioni da effettuare	
            h = waitbar(0,'Please wait...');
            counterr=0;
            for r=1:n2,
                waitbar(r/tot,h)
                s2=0;
                % pongo in s2 la r-sima riga di mat2
                trovazero=find(mat2(r,:)==0);
                if isempty(trovazero),
                    jjj=M2;
                else
                    jjj=trovazero(1,1);								%jjj=campione finale del ciclo
                end;
                s2=mat2(r,1:jjj);
                clear jjj;
                % preparazione del vettore battito tramite sottrazione della retta
                % congiungente i punti estremi (per hp. dipendente dal resp)
                % e normalizzazione tra 0 e 1 delle ampiezze
                J=max(size(s2));                           
                if NORMA==1,
                    % creo la retta che congiunge i punti estremi (primo,ultimo)
                    coangy=(s2(J)-s2(1))/(J-1);	
                    quotay=(s2(1)*J-s2(J))/(J-1);
                    j=linspace(1,J,J);
                    ly=coangy.*j+quotay;                       
                    % sottraggo la retta così che i due punti estremi abbiano
                    % la stessa ampiezza
                    s2=(s2-ly)+max(s2);
                end;
                % derivata prima
                % se devo usarla creo il vettore dy(j)=y(j)-y(j-1)
                y=s2;        
                if (cif==0)|(cif==11)|(cif==1)|(cif==4),
                    dy=zeros(1,J);
                    diffy=diff(y);
                    dy(1,2:J)=diffy;
                    % normalizzo la derivata prima tra 0 e 1
                    dy=(dy-min(dy(2:J)))/(max(dy)-min(dy(2:J)));
                end;
                % se devo usarla, creo la derivata seconda
                if (cif==1)|(cif==2)|(cif==3)|(cif==4),
                    % creo il vettore dy2(j)=y(j)-2*y(j-1)+y(j-2)
                    dy2=zeros(1,J);
                    diffy1=diffy(1,2:length(diffy));
                    dy2(1,3:J)=diffy1-diffy(1,1:length(diffy)-1);   
                    % normalizzo la derivata seconda tra 0 e 1
                    dy2=(dy2-min(dy2(3:J)))/(max(dy2)-min(dy2(3:J)));
                    clear diffy1,clear diffy
                end;
                % normalizzo le ampiezze così che siano tra 0 e 1 
                y=(y-min(y))/(max(y)-min(y));    
                clear ly,clear coangy,clear quotay
                
                %P1NPROVA 
                %P1MOD+DISS
                % Calcolo la matrice d delle dissimilarità locali dei vettori x,y
                % per gli elementi all'interno della finestra di ampiezza r
                % cif=modalità di calcolo della cifra di merito:
                % 0=amp+dx; 10=amp; 11=dx; 1=amp+dx+d2x; 2=amp+d2x
                % 3=d2x;     4=dx+d2x;
                % amp=ampiezza finestra, per segnali di lunghezza uguale, in cui il costo è diverso da 100
                % va inserito dall'utente.
                
                if (I>amp)&(J>amp),
                    rw=amp+abs(I-J);
                    % 17/2/99
                    if (I<=rw),
                        rw=I;
                    end
                    if (J<=rw),
                        rw=J;
                    end;
                else
                    rw=min(I,J);
                end;
                % 6/10/98 evito di calcolare d se la cifra di merito non lo richiede
                if (cif==0)|(cif==10)|(cif==1)|(cif==2),
                    d0=zeros(I,J);
                end
                if (cif==0)|(cif==11)|(cif==1)|(cif==4),
                    d1=zeros(I,J);
                end
                if (cif>=1)&(cif<=4),
                    d2=zeros(I,J);
                end
                % qui inizia il calcolo delle matrici
                for i=1:I              
                    for j=1:J                     
                        if (i<j+rw)&(j<i+rw)
                            if (i==1)|(j==1)
                                if (cif==0)|(cif==10)|(cif==1)|(cif==2),
                                    d0(i,j)=abs(x(i)-y(j));
                                end;
                            end;
                            if (i>1)&(j>1),
                                if (cif==0)|(cif==10)|(cif==1)|(cif==2),
                                    d0(i,j)=abs(x(i)-y(j));
                                end;
                                if (cif==0)|(cif==11)|(cif==1)|(cif==4),
                                    d1(i,j)=abs(dx(i)-dy(j));
                                end;
                                if (i>2)&(j>2)&(cif>=1)&(cif<=4),
                                    d2(i,j)=abs(dx2(i)-dy2(j));
                                end;
                            end;
                        end;
                    end;
                end 
                % costruzione matrice distanza totale come somma matrici ampiezza e derivate
                D=zeros(I,J);
                if (cif==0)|(cif==10)|(cif==1)|(cif==2),
                    D=D+d0;
                    clear d0;
                end;
                if (cif==0)|(cif==1)|(cif==4)|(cif==11),
                    D=D+d1;
                    clear d1;
                end;
                if (cif==1)|(cif==2)|(cif==3)|(cif==4),
                    D=D+d2;
                    clear d2
                end
                d=D;
                clear D
                %--------------------------------------------------------------------------
                % P2MOD
                % inizializzazione della matrice g con vincolo sulla finestra
                i=1;
                j=1;
                g=ones(I,J)*1000;
                %for i=1:I
                %   for j=1:J
                %      if (i>=j+r)|(j>=i+r)
                %         g(i,j)=1000;
                %      end; 
                %  end;
                %end;
                %---------------------------------------------------------------------------
                % PICT3
                % inizializzazione
                i=1;
                j=1;
                g(1,1)=2*d(1,1);
                if p==0
                    cont=1;
                    % cornice orizzontale
                    for j=2:rw,
                        g(1,j)=g(1,j-1)+d(1,j);
                    end
                    
                    % cornice verticale
                    for i=2:rw,
                        g(i,1)=g(i-1,1)+d(i,1);
                    end
                    
                    % interno della cornice
                    i=1;j=2;
                    while j<=J,
                        i=i+1;
                        if i>=j+rw
                            j=j+1;
                            if j>rw
                                i=j-rw;
                            else i=1;
                            end;
                        elseif i<=I
                            g(i,j)=min([g(i-1,j)+d(i,j),
                                g(i-1,j-1)+2*d(i,j),
                                g(i,j-1)+d(i,j)]);
                        end; 
                    end;
                    DI(r)=g(I,J)/(I+J);
                elseif p==1
                    cont=2;
                    % cornice orizzontale
                    for j=2:rw,
                        g(1,j)=g(1,j-1)+d(1,j);
                    end
                    
                    % cornice verticale
                    for i=2:rw,
                        g(i,1)=g(i-1,1)+d(i,1);
                    end
                    % 2° cornice verticale
                    j=2;
                    for i=2:rw,
                        g(i,j)=min([g(i-1,j)+d(i,j),
                            g(i-1,j-1)+2*d(i,j),
                            g(i,j-1)+d(i,j)]);
                    end;
                    % 2° cornice orizzontale
                    i=2;
                    for j=3:rw,
                        g(i,j)=min([g(i-1,j)+d(i,j),
                            g(i-1,j-1)+2*d(i,j),
                            g(i,j-1)+d(i,j)]);
                    end
                    
                    % interno della cornice
                    i=2;j=3;
                    while j<=J,
                        i=i+1;
                        if i>=j+rw
                            j=j+1;
                            if j>rw+2
                                i=j-rw;
                            else i=2;
                            end;
                        elseif i<=I
                            g(i,j)=min([g(i-1,j-2)+2*d(i,j-1)+d(i,j),
                                g(i-1,j-1)+2*d(i,j),
                                g(i-2,j-1)+2*d(i-1,j)+d(i,j)]);
                        end; 
                    end;
                    DI(r)=g(I,J)/(I+J)
                end;
                clear d,clear dy;
                if cif==1,
                    clear dy2
                end;
                %---------------------------------------------------------------------------
                % PICT4
                % inizializzazione puntatori
                % la ricerca procede a ritroso da (I,J) a (1,1)
                i=I;
                j=J;
                k=1;
                WFX(r,k)=i;
                WFY(r,k)=j;
                k=k+1;
                while (i>1)&(j>1),
                    gmin=min([g(i-1,j-1),g(i-1,j),g(i,j-1)]);
                    if gmin==g(i-1,j-1),
                        WFX(r,k)=i-1;
                        WFY(r,k)=j-1;
                        % 17/2/99
                        if (i-1)<=(j-rw)|(j-1)<=(i-rw),
                            %warndlg('Error ! Window too narrow!')
                            counterr=counterr+1;
                        end
                        i=i-1;
                        j=j-1;
                    elseif gmin==g(i-1,j),
                        WFX(r,k)=i-1;
                        if (i-1)<=(j-rw),
                            %warndlg('Error ! Window too narrow!')
                            counterr=counterr+1;
                        end;
                        WFY(r,k)=j;
                        i=i-1;
                    elseif gmin==g(i,j-1),
                        WFX(r,k)=i;
                        WFY(r,k)=j-1;
                        if (j-1)<=(i-rw),
                            %warndlg('Error ! Window too narrow!')
                            counterr=counterr+1;
                        end;
                        j=j-1;
                    end;
                    k=k+1;
                end;
                % sono sulla cornice orizzontale ma non nell'origine
                if (i==1)&(j>1),
                    while j>1,
                        WFX(r,k)=i;
                        WFY(r,k)=j-1;
                        j=j-1;
                        k=k+1;
                    end;
                    % sono sulla cornice verticale ma non nell'origine
                elseif (i>1)&(j==1),
                    while i>1,
                        WFX(r,k)=i-1;
                        WFY(r,k)=j;
                        i=i-1;
                        k=k+1;
                    end;
                end;
                
                %________________________________________________
                % SE VOGLIO VISUALIZZARE DEVO TOGLIERE I COMMENTI
                %plot(wfx,wfy),pause
                % traccio i punti omologhi
                %K=max(size(wfx));
                %points=K;
                %figure(2)
                %plot(1:max(size(s1)),s1(1:max(size(s1))),1:max(size(s2)),s2(1:max(size(s2)))),hold on
                % for i=1:points,
                %	T1=[wfx(i) wfy(i)];
                %	T2=[s1(wfx(i)) s2(wfy(i))];
                %	line(T1,T2)
                %	hold on
                %end
                %pause
                %hold off
                %points=K;
                %figure(3)
                %dif1=[0 diff(s1)];
                %dif2=[0 diff(s2)];
                %plot(1:max(size(dif1)),dif1(1:max(size(dif1))),1:max(size(dif2)),dif2(1:max(size(dif2)))),hold on
                %for i=1:points,
                %	T1=[wfx(i) wfy(i)];
                %	T2=[dif1(wfx(i)) dif2(wfy(i))];
                %	line(T1,T2)
                %	hold on
                %end
                %pause
                %hold off
                % ____________________________
                clear g;
                % passo al successivo battito
                clear s2
            end
            counterr
            caricoWF=1;
            % salva matrici WFX e WFY
            ButtonName=questdlg('Save results on .mat file?','SAVE in Matlab format','Yes');
            if strcmp(ButtonName,'Yes')==1,
                [newfile,newpath]=uiputfile('*.mat','Save Warping Function Matrices and distances');
                if newfile==0,
                    warndlg('No File selected !','!! WARNING !!');
                end;
                save('DTWresults.mat');
            end;
            % Figura 2
            % Visualizzazione serie delle distanze globali dal template (DI)
            figure(2)
            plot(DI),ylabel('Total distance'),xlabel('#beats'),title('TW distance from template')        
        elseif carico==1,
            errordlg('No TW average imported in work area!','!!! WARNING !!!');
        else
            errordlg('No file imported in work area!','!!! WARNING !!!');
        end;
        
        %-----------------------------------------------------------------------------------------
        % questo è l'end del case    
    end
    
    
