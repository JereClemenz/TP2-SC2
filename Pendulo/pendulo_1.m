clc; clear all; close all;

%Ítem [3] Calcular un controlador que haga evolucionar al péndulo en el equilibrio inestable, 
% partiendo de una condición inicial nula en el desplazamiento y termine en -10 metros manteniendo 
% la vertical. Determinar el ángulo máximo que puede alejarse de la vertical en t=0 para que el sistema 
% cumpla el objetivo de control.


m=.1;
Fricc=0.1; 
l=1.6;
g=9.8;
M=1.5;

% Matrices del sistema linealizado
A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 Fricc/(l*M) g*(m+M)/(l*M) 0];
B=[0; 1/M; 0; -1/(l*M)];
C=[1 0 0 0]; 

%Diseño con LQR
Q=1*diag([10 .01 10 1]);    R=5;

%Construcción del Hamiltoniano para el cálculo del controlador
Ha=[A -B*inv(R)*B'; -Q -A'];
[n,va]=size(Ha);
[autovectores,autovalores]=eig(Ha);
MX1X2=[];% obtenemos [M PM]

%Extraigo autovalores negativos
for i=1:1:length(autovalores)
    if (real(autovalores(i,i)))<0
        MX1X2=[MX1X2 autovectores(:,i)];
    end
end    

MX1=MX1X2(1:n/2,:);%M
MX2=MX1X2(n/2+1:end,:);%PM
P=real(MX2*inv(MX1));%P=PMinv(M)

%Calculo de K
K=inv(R)*B'*P;

% K2=lqr(A, B, Q, R) con esta funcion calculo el controlador sin ningun
% calculo adicional o complejo

%Referencia distinta de cero
G=-inv(C*inv(A-B*K)*B);

%Simulación del control:
h=10^-4;%paso
tsim=10; %tiempo de simulacion
t=0:h:(tsim-h);

%Referencia
setpoint_distancia=-10; % la distancia de desplazamiento es -10m 
ref_ang=0;%Idealmente simpre esta en equilibrio inestable

%condiciones iniciales
delta(1)=0;        %x1
delta_p(1)=0;      %x2
theta(1)=0.01;        %x3
theta_p(1)=0;      %x4

estados=[delta(1);
        delta_p(1);
        theta(1);
        theta_p(1)];

Xop=[0 0 0 0]';
x=[delta(1) delta_p(1) theta(1) theta_p(1)]';
theta_pp=0;
for i=1:round(tsim/h)
    
    u(i) = -K*estados+setpoint_distancia*G;
    %Variables del sistema lineal
    delta(i)= x(1);
    delta_p(i)= x(2);
    theta(i)= x(3);
    theta_p(i)= x(4);
    
    %Sistema lineal
%     xp=A*x+B*u1(i);
%     x=x+h*xp;
    
    %Sistema no lineal
    delta_pp=(u(i)-Fricc*x(2)-m*l*theta_pp*cos(x(3)-Xop(3))+m*l*sin(x(3)-Xop(3))*x(4)^2)/(M+m);
    theta_pp=(g*sin(x(3)-Xop(3))-delta_pp*cos(x(3)-Xop(3)))/l;
    
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
end    
    
    
subplot(3, 2, 1);
plot(t,delta);
title('desplazamiento');
xlabel('Tiempo (seg.)');
ylabel('distancia');
grid on;

subplot(3, 2, 2);
plot(t,delta_p);
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