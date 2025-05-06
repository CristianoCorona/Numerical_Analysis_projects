%Definizione parametri
Ibar_s=2e-6;
R=1e3;
C=3e-6;
T=20e-3;
%Variabili ausiliarie
v_inf=R*Ibar_s;
tau=R*C;
G=1/R;
%Definizione funzioni
v = @(t) v_inf*(1-exp(-t/tau));
i_R = @(t) G*v(t);
i_C = @(t) Ibar_s*exp(-t/tau);

%1) DERIVAZIONE NUMERICA
%Vettore distanze h sempre più piccole
hvec = T./(2.^[0:10]');
%Vettore errori associati alla distanza h
errvec=zeros(length(hvec),1);

for i=1:length(hvec)
    %Prendo la distanza i-esima
    h=hvec(i);
    %Vettore nodi (istanti temporali)
    tn=[0:h:T]';
    %Tensioni e correnti per ciascun istante
    vn=v(tn);
    iCn=i_C(tn);

    dvdt=zeros(length(tn),1);
    %dvdt:differenza in avanti
    for j=1:length(tn)-1
        dvdt(j)=(vn(j+1)-vn(j))/h;
    end
    %differenza all'indietro
    dvdt(length(tn))=(vn(length(tn))-vn(length(tn)-1))/h;
    %approssimazione iCh
    iCh=C*dvdt;
    %errore nodale per h
    en = iCn - iCh;
    errvec(i)=max(abs(en));
end 
setfonts();
[p,C]=order_estimate(hvec,errvec);
format short e;
