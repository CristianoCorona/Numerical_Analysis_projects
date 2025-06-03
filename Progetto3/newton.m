function [x_vect,it] = newton(x0,nmax,toll,fun,J)
% [x_vect,it] = newton(x0,nmax,toll,fun,J)
%
% Metodo di Newton per la ricerca degli zeri della funzione vettoriale fun.
% Test d'arresto basato sull'incremento tra due iterazioni successive.
%
% Parametri di ingresso:
% x0 Punto di partenza
% nmax Numero massimo di iterazioni
% toll Tolleranza sul test d'arresto
% fun dfun Anonymous function contenenti la funzione e la sua derivata
%
% Parametri di uscita:
% x_vect Vett. contenente tutte le iterate calcolate
% it Iterazioni effettuate

x_vect = x0;
err = toll+1;
it = 0;

while it < nmax && err > toll
    x_old = x_vect;
    fx = fun(x_old);
    Jx = J(x_old);
    if det(Jx) == 0
        % error fa terminare l'esecuzione della function
        error(' Arresto per sistema singolare');
    end
    % nuova iterazione di newton
    x_vect = x_old - Jx\fx;
    err = norm(x_vect-x_old,Inf);
    % salvo le quantita` di interesse negli output
    it = it+1;
end
return
