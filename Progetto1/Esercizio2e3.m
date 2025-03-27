%Esercizio 2:Metodi diretti
%1) matrice ammettenze Y e termine noto b
Vin=5;
Rin=600;
R1=Rin/1;
R2=Rin/2;
R3=Rin/3;
R4=Rin/4;

Gin=1/Rin;
G1=1/R1;
G2=1/R2;
G3=1/R3;
G4=1/R4;

n=6;
%starting struct
Y=zeros(n)+G2*diag(ones(n-3,1),3)+G2*diag(ones(n-3,1),-3);
%details
Y(1,2) = G1;
Y(2,1) = G1;
Y(2,3) = G1;
Y(3,2)=G1;
Y(4,5)=G3;
Y(5,4)=G3;
Y(5,6)=G3;
Y(6,5)=G3;
diagonal=[-(Gin+G1+G2),-(2*G1+G2),-(G1+G2),-(G2+G3+G4),-(G2+2*G3+G4),-(G2+G3+G4)];
Y=Y+diag(diagonal);
%la matrice sarebbe così definita negativa; ha senso moltiplicare tutte le
%equazioni per -1, ottenendo una matrice SDP. RICORDA ANCHE b!
Y=-Y;
b=zeros(n,1);
b(1,1)=Gin*Vin; %da modficare qualora non si cambiasse il segno!

format short e;
%2) Verifica esistenza e unicità fattorizzazione LU
%Sfrutto la condizione sufficiente; Y è SDP (simmetrica definita positiva)
if(Y == Y')
    disp('matrice simmetrica!');
else
    disp('matrice non simmetrica :(');
end

eigen=eig(Y);
for i=eigen
    if(i<=0)
        disp('matrice non SDP');
    else 
        disp('matrice SDP! Esiste la fattorizzazione LU!');
    end
end

%3)Calcolo fattorizzazione lu
[L,U,P] = lu(Y);

%4) Risoluzione sistema
format long e;
y = fwsub(L,P*b); %in questo caso P non ha effetto, ma cambiando i parametri del problema potrei avere pivoting!
xc = bksub(U,y);
xm = Y\b;
%5)
format short e;
err = norm(xm-xc,inf);

%Esercizio 3: metodi iterativi
%1) Verifica Y SDP; già verificata al punto 2.2
%2) Soluzione con GS
format long e;
x0 = zeros(6,1);
toll = 1e-12;
nmax = 1000;

[xGS,kGS] = gs(Y,b,x0,toll,nmax);

%3) errore relativo GS
format short e;
true_rel_err_GS = norm(xm-xGS)/norm(xm);

%4)raggio spettrale e confronto con errore relativo
%so già dal punto precedente che la matrice rispetta le ipotesi di
%applicabilità del metodo

D = diag(diag(Y));
E = -tril(Y,-1);
F = -triu(Y,1);
Bgs=(D-E)\F;
rhoGS = max(abs(eig(Bgs)));
rhoGSallaK=rhoGS^kGS;

%5) Richardson
format long e;
alpha_opt = 2/(min(eig(Y))+max(eig(Y)));
[xR,kR] = richardson(Y,b,x0,alpha_opt,toll,nmax);

true_err_rel_R = norm(xm-xR)/norm(xm);
Br = eye(n)-alpha_opt*Y;
rhoR = max(abs(eig(Br)));
rhoRallak = rhoR^kR;




