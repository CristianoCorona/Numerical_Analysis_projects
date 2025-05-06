function [xvect,it] = newton (x0, nmax, toll, fun, dfun)

%Metodo di Newton per la ricerca di radici della funzione fun

%x0: punto di partenza
%nmax: numero massimo di iterazioni
%toll: tolleranza sul test d'arresto
%fun: funzione (function handle)
%dfun: derivata delal funzione (function handle)

%xvect: vettore contenente le soluzioni a ogni iterazione (l'ultima componente sarà la soluzione
%it : numero di iterazioni effettuate

%Inizializzazione
it = 0;
incr = toll+1; %inizializzata per entrare nel while
xvect=x0;

while(it < nmax && incr >=toll)
    xv = xvect(end);
    if(dfun(xv)==0)
        disp ('Arresto per annullamento dfun')
        break
    else
        xn = xv - fun(xv)/dfun(xv);
        incr = abs(xn-xv);
        xvect = [xvect;xn];
        it = it +1;
    end
end
fprintf('Numero di iterazioni: %d\n',it);
fprintf('Zero calcolato: %12.13f \n',xvect(end));
return
