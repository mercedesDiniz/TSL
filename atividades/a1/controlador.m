%% Atividade A1: Estudo de caso para o controle de velocidade longitudinal de um quadricóptero.
clear all; close all; clc;

%% Setup
% Carregando os coeficientes do modelo de 2ª ordem
load('sys_mod01.mat');
Ts = model.Ts;
Az = model.Az; a1 = Az(2); a2 = Az(3);
Bz = model.Bz; b0 = Bz(2); b1 = Bz(3);

% Difinindo a função de transferencia da planta
Gz = tf(Bz,Az,Ts);

%% Controlador PID
% Utilizando o modelo de 2ª ordem para sintonia
kp = 0.2; ki = 2; kd = 0.1;

% PID digital via aproximação Backward Difference
s0 = kp +ki*Ts +kd/Ts;
s1 = -kp -2*kd/Ts;
s2 = kd/Ts;

% Difinindo a função de transferencia do controlador
% C(z) = s0 + s1*z^(-1) + s2*z^(-2) /  1 + z^(-1)
Cz1 = tf([s0 s1 s2],[1  -1  0],Ts);

% Salvando o modelo
control.Cz = Cz1;
control.Ts = Ts;
save 'control_pid.mat' control;

%% Controlador PI-D
% Mesmos ganhos kp e ki do controlador PID
kd = 0; % resultando em uma estrutura PI

% PID digital via aproximação Backward Difference
s0 = kp +ki*Ts +kd/Ts;
s1 = -kp -2*kd/Ts;
s2 = kd/Ts;

% Difinindo a função de transferencia do controlador
% C(z) = s0 + s1*z^(-1) + s2*z^(-2) /  1 + z^(-1)
Cz2 = tf([s0 s1 s2],[1  -1  0],Ts);

% Salvando o modelo
control.Cz = Cz2;
control.Ts = Ts;
save 'control_pi-d.mat' control;

