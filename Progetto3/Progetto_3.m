%%Progetto 3

%%%%
%Esercizio 2: Modello Numerico
constants;
%altre costanti ausiliarie
R = 20;
i_s = 0.1;
i_sat = 1e-9;
A = (2/3)*R*i_sat;
B = A + (R*i_s)/3;

%1) Grafico
x_a  = 0.4;
x_b = 0.5;

xx = [x_a:(x_b-x_a)/10000:x_b];
g_1 = @(x) x;
g_2 = @(x) B-A*exp(x./Vth);

g_1dis = g_1(xx);
g_2dis = g_2(xx);

figure(1)
plot(xx,g_1dis,"b",xx,g_2dis,"r","LineWidth",2);
title('g_1(x) e g_2(x)');
xlabel("x")
grid on;
legend("g_1(x)","g_2(x)");
exportgraphics(gcf, 'GraficoG1eG2.pdf', 'ContentType', 'vector', ...
               'BackgroundColor', 'white', 'Resolution', 1200);

%2) x* con fzero
x0=0;
g = @(x) g_1(x)-g_2(x);
xex = fzero(g,x0);
format long e;
fprintf("xex = %d\n",xex);

%3) Iterazioni di punto fisso: previsioni teoriche
Tg = @(x) Vth*log((B-x)./A);
dTg = @(x) -Vth./(B-x);
W = abs(dTg(xex));
fprintf("W = %d\n",W);

%4) Iterazioni punto fisso
x0 = 0;
itmax = 1000;
toll = 1e-12;
[x, niter, err] = fixed_point(x0, Tg, toll, itmax);
x_fix = x(end);
EST_ERR = err(end);
TRUE_ERR = xex - x_fix;

fprintf("x_fix = %d\n",x_fix);
fprintf("EST_ERR = %d\n",EST_ERR);
fprintf("TRUE_ERR = %d\n",TRUE_ERR);

%3) Regime Dinamico
C = 1e-6;
is_bar =  0.1;
freq = 100;
omega = 2*pi*freq;
td = 1e-2;

t0 = 0;
tf = 0.1;
tspan = [t0 tf];
y0 = 0;

is_fun = @(t) is_bar.*sin(omega.*t).*exp(-t./td);

%A) ode15s
f = @(t,x) -(3/(2*R*C)).*x + (i_sat/C).*(1-exp(x./Vth))+is_fun(t)/(2*C);
[tm, xm] = ode15s(f,tspan,y0);
delta_xA = max(xm) - min(xm);
%B) crank_nicolson
dfdx = @(t,x) -3/(2*R*C) - (i_sat/(C*Vth))*exp(x./Vth);
NT = 1000;
h = (tf-t0)/NT;
tol = 1e-12; %sono hardcoded nella funzione crank_nicolson, non vengono
% % % % % % maxit = 100; %presi come input
[tth, xth] = crank_nicolson(f,dfdx,t0,tf,y0,h);
delta_xB = max(xth) - min (xth);

fprintf("delta_xA = %d\n",delta_xA);
fprintf("delta_xB = %d\n",delta_xB);

figure(2)
plot(tm,xm,"b-",tth,xth,"r--","LineWidth",2);
title('Soluzioni ODE: ode15s vs Crank-Nicolson');
xlabel('Tempo (s)');
ylabel('x(t)');
grid on;
legend("ode15s","Crank-Nicolson");
exportgraphics(gcf, 'GraficoODE.pdf', 'ContentType', 'vector', ...
               'BackgroundColor', 'white', 'Resolution', 1200);