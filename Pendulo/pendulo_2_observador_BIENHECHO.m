clc; clear all; close all;

% Ítem [4] Incorporar un observador para el caso en que sólo puedan medirse el desplazamiento delta
% y el ángulo phi, repetir las simulaciones para las condiciones anteriores y graficar los resultados en 
% gráficas superpuestas.

%Version 2. Cambios sugeridos por Pucheta --> C con la fila de salida de
%angulo

m=.1;
Fricc=0.1; 
l=1.6;
g=9.8;
M=1.5;

% Matrices del sistema linealizado
A=[0 1 0 0;0 -Fricc/M -m*g/M 0; 0 0 0 1; 0 Fricc/(l*M) g*(m+M)/(l*M) 0];
B=[0; 1/M; 0; -1/(l*M)];
C=[1 0 0 0;0 0 1 0]; 

%Diseño con LQR
Q=1*diag([10 .01 10 1]);    R=5;
K=lqr(A, B, Q, R)

%Diseño con LQR para el observador
% Qo=1*diag([1 1 1 1]);    Ro=diag([1 1]);
% Qo=1e4*diag([1 10 1 10]);    Ro=0.001
Qo=1*diag([1 10 1 10]);    Ro=1;
Ao=A';
Bo=C';
Co=B';
%controlador del observador calculado con funcion 'lqr'
Ko=lqr(Ao,Bo, Qo, Ro);

%Referencia distinta de cero
G=-inv(C(1,:)*inv(A-B*K)*B);

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
theta(1)=0.1;     %x3
theta_p(1)=0;      %x4

estados=[delta(1);
        delta_p(1);
        theta(1);
        theta_p(1)];
    
estados_obs=[delta(1);
        delta_p(1);
        theta(1);
        theta_p(1)];
    
Xop=[0 0 0 0]';
x=[delta(1) delta_p(1) theta(1) theta_p(1)]';
%inicializacion del observador
xobs=[0 0 0 0]';

theta_pp=0;

for i=1:round(tsim/h)
    u(i)= -K*estados_obs+G*setpoint_distancia;

    %Variables del sistema
    delta(i)= x(1);
    delta_p(i)= x(2);
    theta(i)= x(3);
    theta_p(i)= x(4);
    
    %Sistema no lineal
    delta_pp=(u(i)-Fricc*x(2)-m*l*theta_pp*cos(x(3))+m*l*sin(x(3))*x(4)^2)/(M+m);
    theta_pp=(g*sin(x(3))-delta_pp*cos(x(3)))/l;
    
    x_p_1=x(2);
    x_p_2=delta_pp;
    x_p_3=x(4);
    x_p_4=theta_pp;
    
    xp=[x_p_1;x_p_2;x_p_3;x_p_4];
%     xp=A*x+B*u(i);
    x=x+h*xp;
    
    %------con Observador----------------------
  
    delta_o(i)= xobs(1);
    delta_p_o(i)= xobs(2);
    theta_o(i)= xobs(3);
    theta_p_o(i)= xobs(4);
    
    y_sal_o = C * estados_obs;
    y_sal   = C * estados;
    
    e=y_sal-y_sal_o;
    
    x_antp     = A*xobs+B*u(i)+Ko'*e;%complicacion
    xobs       = xobs + x_antp*h;
    
    %-------------------------------------------
    
    estados=[delta(i);
        delta_p(i);
        theta(i);
        theta_p(i)];
    
    estados_obs=[delta_o(i);
        delta_p_o(i);
        theta_o(i);
        theta_p_o(i)];
    
end    
    
    
%_______________SIN OBS---------------------------------------------------
delta_so(1)=0;        %x1
delta_p_so(1)=0;      %x2
theta_so(1)=theta(1);     %x3
theta_p_so(1)=0;      %x4
estados_so=[delta_so(1);
        delta_p_so(1);
        theta_so(1);
        theta_p_so(1)];

Xop_so=[0 0 0 0]';
x_so=[delta_so(1) delta_p_so(1) theta_so(1) theta_p_so(1)]';
theta_pp_so=0;
for i=1:round(tsim/h)
    
    %---------SIN OBS---------------------------
    u1_so(i) = -K*estados_so+setpoint_distancia*G;
    %Variables del sistema lineal
    delta_so(i)= x_so(1);
    delta_p_so(i)= x_so(2);
    theta_so(i)= x_so(3);
    theta_p_so(i)= x_so(4);
    
    %Sistema lineal sin obs
%     xp_so=A*x_so+B*u1_so(i);
    delta_pp_so=(u1_so(i)-Fricc*x_so(2)-m*l*theta_pp_so*cos(x_so(3))+m*l*sin(x_so(3))*x_so(4)^2)/(M+m);
    theta_pp_so=(g*sin(x_so(3))-delta_pp_so*cos(x_so(3)))/l;
    
    x_p_1_so=x_so(2);
    x_p_2_so=delta_pp_so;
    x_p_3_so=x_so(4);
    x_p_4_so=theta_pp_so;
    
    xp_so=[x_p_1_so;x_p_2_so;x_p_3_so;x_p_4_so];
    x_so=x_so+h*xp_so;
    
    estados_so=[delta_so(i);
        delta_p_so(i);
        theta_so(i);
        theta_p_so(i)];
    
end    
    
%----------------------------------------------------------------------


subplot(3, 2, 1);
hold on
plot(t,delta);
plot(t,delta_so);
hold off
legend({'Con observador','Sin observador'})
title('desplazamiento');
xlabel('Tiempo (seg.)');
ylabel('distancia');
grid on;

subplot(3, 2, 2);
hold on
plot(t,delta_p);
plot(t,delta_p_so);
hold off
legend({'Con observador','Sin observador'})
title('Velocidad');
xlabel('Tiempo (seg.)');
ylabel('Velocidad (m/s)');
grid on;

subplot(3, 2, 3);
hold on
plot(t,theta*(180/pi));
plot(t,theta_so*(180/pi));
hold off
legend({'Con observador','Sin observador'})
title('Poscion angular \theta_t');
xlabel('Tiempo (seg.)');
ylabel('Posicion angular (Rad)');
grid on;

subplot(3, 2, 4);
hold on
plot(t,theta_p);
plot(t,theta_p_so);
hold off
legend({'Con observador','Sin observador'})
title('Velocidad angular \omega_t');
xlabel('Tiempo (seg.)');
ylabel('Posicion angular (Rad)');
grid on;

subplot(3, 1, 3);
hold on
plot(t,u);
plot(t,u1_so);
hold off
legend({'Con observador','Sin observador'})
title('Accion de control u_t');
xlabel('Tiempo (seg.)');
ylabel('V');
grid on;