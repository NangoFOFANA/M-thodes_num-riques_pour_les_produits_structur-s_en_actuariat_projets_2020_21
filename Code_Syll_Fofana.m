%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% CODE projet Méthode numériques pour les produits structurés en
%%%% actuariat SYLL & FOFANA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%  parametre généraux   %%%%%%%%%%%

%Parametres generaux
clear all
close all
nx=21;
ny=21;
rhoxr=0;
R1=0.05;
n=(nx-1)*(ny-1);
T=10;
eta=2;
delta=0.03;
TMG=0.025;
sigmax=0.008;
sigmar=0.02;
rhoax=0.95;
rhoar=0.25;
mu_a=0.04;
sigmaA=0.06;
kx=0.3;
kr=0.18;
rho=0.95;
r0=0.03;
rinf=0.05;
x0=0.005;
xinf=0.005;
mu_i=0.04;
L1=3*sigmax/(sqrt(2*kx));
xmin=xinf-L1;
xmax=xinf+L1;
Dx=(xmax-xmin)/nx;
xx = linspace(xmin,xmax,n);
L2=3*sigmar/(sqrt(2*kr));
rmin=rinf-L2;
rmax=rinf+L2;
Dr=(rmax-rmin)/ny;
rr = linspace(rmin,rmax,n);
T=10;
N=21;
Dt=T/N;
tt = linspace(0,T,N);
%% Resolution de notre EDP phi(n+1)=inv(A)*(1/Dt*phi(n)+g(x))
%%%%% système de l'EDP
%%%%% definition des parametres du vecteur phi
D=-rho*sigmax*sigmar/(Dx*Dr);
B=-(sigmar)^2/(2*(Dr)^2);
C=-(sigmax)^2/(2*(Dx)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Matrice(A) du system de l'edp
A=zeros(n);
for i=nx:n
    A(i,i-(nx-1))=B;
end
%%%%%%%%%%%%%
for i=2:n
    A(i,i-1)=C;
end
%%%%%%%%%%%%%
for i=1:((nx-2)*(ny-1)-1)
    A(i,i+nx)=D;
end
%%%%%%%
for i=1:(n-1)
    if (mod(i,nx)==0)
        A(i,i+1)=-(kx/Dx*(xinf-xx(i))-B+D);
        A(i+1,i)=(kx/Dx*(xinf-xx(i))-B+D);
    end
end
%%%%%
for i=1:((nx-2)*(ny-1))
    A(i,i+(nx-1))=-((kr/Dr)*(rinf-rr(i))-C+D);
    A(i+(nx-1),i)=((kr/Dr)*(rinf-rr(i))-C+D);
    
end
%%%%%%%
for i=1:n
    if ((mod(i,nx)==0))
        A(i,i)=1/Dt+(kx/Dx)*(xinf-xx(i))+(kr/Dr)*(rinf-rr(i))-2*B-2*C+D-(gonc(xx(i),rr(i))-rr(i)-fonc(xx(i)));
    end
end

A(1,1)=1;
A(1,2)=2; 
A(1,3)=3;
A(1,ny)=0;
A(nx,1)=0;
A(n,n)=1/Dt+(kx/Dx)*(xinf-xx(n))+(kr/Dr)*(rinf-rr(n))-2*B-2*C+D-(gonc(xx(i),rr(n))-rr(n)-fonc(xx(n)));
%%%% resolution de du système de l'EDP

phi=ones(n,1);

%D=1/Dt*U+g;
%b=zeros(n,1)
for t=1:N
    b=1/Dt*phi+fonc(xx)';
    phi=A\b;%inv(A)*b';
end

%%% Graphique de phi en 2D
PHI=ones(n,1);
for t=1:N
    J=1/Dt*PHI+fonc(xx);
    PHI=A\J;%inv(A)*b';
end
LLL=surf(PHI);
legend([LLL],'Graphique de phi');
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Simulation  par monte carlo
M=10^4; %Nombres de simulation
Z=randn(M); %%% 
%%%% R1
for p=1:M
    r1(p)=r0*exp(-kr)+rinf*(1-exp(-kr))+sqrt(sigmar^2*((1-exp(-2*kr))/2*kr))*Z(p);
end

%%%%%%%%%%%%
%% loi normale centré reduite
W=randn(M); 
V=randn(M);
%%%%%% Fonction x1
for p=1:M
    x1(p)=x0*exp(-kx)+xinf*(1-exp(-kx))+sqrt(sigmax^2*((1-exp(-2*kx))/2*kx))*(rhoax/(sqrt(1-rhoar^2))*V(p)+sqrt((1-rhoar^2-rhoxr^2)/(1-rhoar^2))*W(p));

%x0+kx*(xinf-x0)+(rhoax*sigmax/(sqrt(1+(rhoar)^2)))*W(p)+(sqrt(1-(rhoax)^2-(rhoar)^2/(1-(rhoax)^2)))*V(p);
end
%%%%% fonction RA
ZA=randn(M);
for p=1:M
    rA(p)=mu_a+sigmaA*(rhoar*Z(p)+sqrt(1-rhoar^2)*ZA(p));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LES passifs
%%Fonds propres à la date t=1
PM0=36000000;
F1=PM0*fonc(x0)*Dt;

%% Provision mathématiques à la date t=1
PM1=PM0*exp(Dt*(gonc(x0,r0)-fonc(x0))); %%% fonction PM1

%%les actifs
A0=10147920;
for p=1:M
    A1(p)=A0*(1+rA(p));
end
%%%%% histogramme des Actifs
K=histogram(A1)
legend([K],'Histogramme de A1')
%%%%%%%%%%%%%%%%
L1=phi(2);
L0=phi(1);
%%%%%% Best estimate à la date t=1
%BE1=zeros(M,1);
for t=1:N-1
    BE(t)=PM1*phi(t);%*(x1);
end

BE1=PM1*phi(N-1);
%%% moyenne du best Estimate
moyenne=mean(BE);
%%% Ecart type du best Estimate
sigma=std(BE);
   
%%%  Dette vis à vis des assurés
for p=1:M
    E1(p)= A1(p)-F1-BE1;
end


%%%% Calcul du SCR0
T=ceil((M/200));
VarE1=E1(T);
E0=A0-L0;
SCR0=E0-exp(-r0*Dt)*VarE1;

%%% histogramme de E1
h=bar(E1)
legend([h],'Histogramme de E1')
%%%%%
%% Taux de couverture
TAUCOUV=E0/SCR0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%ù