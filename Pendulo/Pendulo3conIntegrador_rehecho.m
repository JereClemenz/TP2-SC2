clc; clear all; close all;

%Este es el punto 3 correctamente implementado. No resto al modelo liuneal
%y controladores calculados con modelos


m=.1;
Fricc=0.1; 
l=1.6;
g=9.8;
M=1.5;

% Matrices del sistema linealizado equilibrio estable
A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 -Fricc/(l*M) -g*(m+M)/(l*M) 0];
B=[0; 1/M; 0; 1/(l*M)];
C=[1 0 0 0]; 
%usar modelo linealizado 180/pi


% Construcción del sistema ampliado
A=[A zeros(4,1); -C 0];
B=[B; 0];
Ca=[C 0];

%Diseño con LQR
% Q=1*diag([10 .01 10 1]);    R=5;
Q=1*diag([0.1 100 10 10000 1]);    R=100;
Kamp=lqr(A,B,Q,R); %controlador para primera etapa hasta 2 metros sin carga
%ampliado por el integrador
K1=Kamp(1:4);
K1I=-Kamp(5);


%Sistema con 10 veces mas masa---------------------------------------------
m=.1*10;
Fricc=0.1; 
l=1.6;
g=9.8;
M=1.5;
Am=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 -Fricc/(l*M) -g*(m+M)/(l*M) 0];
Bm=[0; 1/M; 0; 1/(l*M)];
Cm=[1 0 0 0]; 

% Construcción del sistema ampliado
Am=[Am zeros(4,1); -Cm 0];
Bm=[Bm; 0];
% C=[C 0];

%Diseño con LQR
% Q=1*diag([10 .01 10 1]);    R=5;
Qm=1*diag([1 1000 1 10000 1]);    Rm=100;
Kamp_m=lqr(Am,Bm,Qm,Rm); %controlador para primera etapa hasta 2 metros sin carga
%ampliado por el integrador
K1m=Kamp_m(1:4);
K1Im=-Kamp_m(5);
%--------------------------------------------------------------------------



%Simulación del control:
h=10e-5;%paso
tsim=100; %tiempo de simulacion
t=0:h:(tsim-h);

%Referencia
% sp_dist=2; % la distancia de desplazamiento es -10m 
%Nueva referencia que en 10 segundos va a 2m y despues vuelve a 0m en 10s
sp_dist=1*square(2*pi*t/tsim)+1;
ref_ang=pi;
pasos=round(tsim/h);
m=ones(1,pasos);
m=m*0.1;
m((pasos/2):end)=m((pasos/2):end)*10;


%condiciones iniciales
delta(1)=0;        %x1
delta_p(1)=0;      %x2
theta(1)=pi;     %x3
theta_p(1)=0;      %x4
psi(1)=0;

estados=[delta(1);
        delta_p(1);
        theta(1);
        theta_p(1)];
integracion(1)=psi(1);    
Xop=[0 0 pi 0]';
x=[delta(1) delta_p(1) theta(1) theta_p(1)]';

theta_pp=0;

%Inicializo los controladores
K=K1;KI=K1I;
for i=1:round(tsim/h)
    
    if m(i)>=0.5
        K=K1m;
        KI=K1Im;
    end
    
    psi_p=sp_dist(i)-C*(estados-Xop);
    psi(i)=integracion+psi_p*h;
    
    u(i)= -K*(estados(1:4)-Xop(1:4))+KI*psi(i);
    
    %Variables del sistema lineal
    delta(i)= x(1);
    delta_p(i)= x(2);
    theta(i)= x(3);
    theta_p(i)= x(4);
    
    
    delta_pp=(u(i)-Fricc*x(2)-m(i)*l*theta_pp*cos(x(3))+m(i)*l*sin(x(3))*x(4)^2)/(M+m(i));
    theta_pp=((g*sin(x(3)))-delta_pp*cos(x(3)))/l;%no restarle 
    x_p_1=x(2);
    x_p_2=delta_pp;
    x_p_3=x(4);
    x_p_4=theta_pp;
    xp=[x_p_1;x_p_2;x_p_3;x_p_4];
    x=x+h*xp;
    
    estados=[delta(i);
        delta_p(i);
        theta(i);
        theta_p(i)];
    integracion=psi(i);
    
end
    
    
subplot(3, 2, 1);
hold on
plot(t,delta);
plot(t,sp_dist,'g--');
hold off
title('desplazamiento');
xlabel('Tiempo (seg.)');
ylabel('distancia');
grid on;

subplot(3, 2, 2);
hold on
plot(t,delta_p);
hold off
title('Velocidad');
xlabel('Tiempo (seg.)');
ylabel('Velocidad (m/s)');
grid on;

subplot(3, 2, 3);
hold on
plot(t,theta*(180/pi));
hold off
title('Poscion angular \theta_t');
xlabel('Tiempo (seg.)');
ylabel('Posicion angular (Rad)');
grid on;

subplot(3, 2, 4);
hold on
plot(t,theta_p);
hold off
title('Velocidad angular \omega_t');
xlabel('Tiempo (seg.)');
ylabel('Posicion angular (Rad)');
grid on;

subplot(3, 1, 3);
hold on
plot(t,u);
hold off
title('Accion de control u_t');
xlabel('Tiempo (seg.)');
ylabel('V');
grid on;
