function [t_h,U,iter_nwt]=crank_nicolson(f,df,t0,T,y_0,h)

%
% Risolve il problema di Cauchy
%
% y'=f(y,t)
% y(0)=y_0
%
% utilizzando il metodo di Crank-Nicolson.
% Per ciascun istante temporale si calcola u_(n+1)
% trovando lo zero dell'equazione (a priori non lineare)
% F(u) = u - u_n - (h/2)*[f(t_n,u_n) + f(t_{n+1},u)] = 0
% con il metodo di Newton. Per questo è necessario passare
% in input l'espressione di df/dy, in modo da creare
% dF(u)/du = 1 - (h/2)* df/dy(t_{n+1},u).
%
% Input:
% -> f: function che descrive il problema di Cauchy
%       (dichiarata come anonymous function)
%       deve ricevere in ingresso due argomenti: f = @(t,y)
% -> df: function che descrive df/dy, anch'essa anonymous function che riceve in
%        ingresso 2 argomenti: df = @(t,y)
% -> t0: istante iniziale
% -> T: l'istante finale dell' intervallo temporale di soluzione
% -> y_0: il dato iniziale del problema di Cauchy
% -> h: l'ampiezza del passo di discretizzazione temporale.
%
% Output:
% -> t_h: vettore degli istanti in cui si calcola la soluzione discreta
% -> U: la soluzione discreta calcolata nei nodi temporali t_h
% -> iter_nwt: vettore che contiene il numero di iterazioni
%                del metodo di Newton necessarie a risolvere l'equazione
%                non lineare ad ogni istante temporale.

t_h = t0:h:T;
% inizializzo il vettore che conterra' la soluzione discreta
N   = length(t_h);
U = zeros(1,N);
U(1) = y_0;
% parametri per il metodo di Newton
nmax = 100;
toll = 1e-12;
iter_nwt=zeros(1,N);
for it=2:N
    
    t_old = t_h(it-1);
    u_old = U(it-1);
    t_nwt = t_h(it);
    
    % Funzioni per il metodo di Newton:
    phi  = @(u) u - u_old - (h/2)*( f(t_old,u_old) + f(t_nwt, u) );
    dphi = @(u) 1 - (h/2)*df(t_nwt, u);
    
    % sottoiterazioni
    [u_nwt, it_nwt] = newton(u_old, nmax, toll, phi, dphi);
    U(it) = u_nwt(end);
    
    % tengo traccia dei valori di it_nwt per valutare la convergenza 
    % delle iterazioni di newton
    iter_nwt(it) = it_nwt;
    
end
