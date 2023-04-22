clc;clear all;close all;


%Parametros de motor
Laa=366e-6;J=5e-9;
Ra=55.6;Bm=0;
Ki=6.49e-3;
Km=6.53e-3;


%Matrices de estado
A = [-Ra/Laa -Km/Laa 0;Ki/J -Bm/J 0;0 1 0];
B = [1/Laa 0;0 -1/J;0 0]; %considerando el torque
C = [0 0 1];             %salida posicion
D = [0 0];


%------Controlador--------------
M = [B(:,1) A*B(:,1) A^2*B(:,1)]; %matriz de controlabilidad solo de Va
rank(M)                 %chequea el rango de la matris M, tiene que ser 3

% Calculo de W    
polAcar = poly(A);  %polinomio caracteristico de A -> |sI-A|
W = [polAcar(3) polAcar(2) 1;
     polAcar(2)      1     0;
         1           0     0];   
   
%Transformacion matricial
T = M * W;

A_controlable = inv(T) * A * T; 
B_controlable = inv(T) * B(:,1);


%  ubicacion de los polos a lazo cerrado
disp('Polos a lazo abierto: ')
eig(A)%polos a lazo abierto, estos son los polos que puedo mover



%seleccion de polos de lazo cerrado (5-26) alfas que me permiten encontrar
%deltas
p1=-500;p2=-1000;p3=-3.0e3;
alfa_i = poly([p1 p2 p3]);%valores caracteristicos del polinomio a buscar ALFAS

%calculo del controlador 
K = (fliplr(alfa_i(2:end)-polAcar(2:end))*inv(T));


disp('Polos a lazo cerrado: ')
eig(A-B(:,1)*K)

%agregamos el efecto del torque formula (5-54)/clase 5
G =-inv(C*inv(A-B(:,1)*K)*B(:,1));

%Variables
tsim=10; 
h=1e-6; 
t=0:h:(tsim-h);

%Entradas del sistema
flag=1;
contador=0;
ref=zeros(1,round(tsim/h));
for i=(round(1/h)):1:(tsim/h)
    if flag==1 
    ref(1,i)=pi/2;
    contador=contador+1;
        if contador==round(2/h)
            flag=0;
            contador=0;
        end
    else 
    ref(1,i)=-pi/2;
    contador=contador+1;
             if contador==round(2/h)
                flag=1;
                contador=0;
            end
    end
   
end
% figure(1)
% plot(t,ref);
% title('Referencia');


flag=1;
contador=0;
tLin=zeros(1,round(tsim/h));
for i=(round(1/h)):1:(tsim/h)
    if flag==1 
    tLin(1,i)=1.15e-3*.1;
    contador=contador+1;
        if contador==round(2/h)
            flag=0;
            contador=0;
        end
    else 
    tLin(1,i)=-0;
    contador=contador+1;
             if contador==round(2/h)
                flag=1;
                contador=0;
            end
    end
   
end
% figure(1)
% plot(t,tLin);
% title('Torque de entrada');

%condiciones iniciales
ia(1)=0;                %Corriente de armadura x1
theta(1)=0;             %Valocidad angular x2
omega(1)=0;             %Posicion angular x3

estados=[ia(1);
        omega(1);
        theta(1)];

    
Xop=[0 0 0]';
x=[ia(1) omega(1) theta(1)]';


B1=[1/Laa;0;0];
B2=[0;-1/J;0];
for i=1:1:(tsim/h)
    
    u1(i) = -K*estados+G*ref(i);
    %Variables del sistema lineal
    ia(i)= x(1,1);
    theta(i)= x(2,1);
    omega(i)= x(3,1);
    
    %Sistema lineal
    xp=A*x+B1*u1(1,i)+B2*tLin(i);
    x=x+h*xp;
    estados=[ia(i);omega(i);theta(i)];
end


% figure(3)
% plot(t,ia);
% title('Corriente de armadura i_a');
% xlabel('Tiempo (seg.)');
% ylabel('Corriente (A)');
% grid on;
% 
% figure(4)
% plot(t,omega);
% title('Velocidad angular \omega_r');
% xlabel('Tiempo (seg.)');
% ylabel('Velocidad angular (rad/s)');
% grid on;

figure(5)
hold on
plot(t,theta);
plot(t,ref);
hold off
title('Poscion angular \theta_t');
xlabel('Tiempo (seg.)');
ylabel('Posicion angular (Rad)');
grid on;
