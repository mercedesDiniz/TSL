%% Atividade 02: Controle do pêndulo invertido sobre carro usando LQG.
clear all; close all; clc;

%% Planta do pêndulo invertido 
% Ref.: https://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling 
% A variável sensorialmente medida é a posição angular do pêndulo (phi)

% Parâmetros do sistema
M = 0.5; m = 0.2; b = 0.1; I = 0.006; g = 9.8; l = 0.3;   

% Modelo em Espaço de Estados Contínuo
p = I*(M+m)+M*m*l^2; 

A = [0      1              0           0;
     0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
     0      0              0           1;
     0 -(m*l*b)/p       m*g*l*(M+m)/p  0];

B = [0; (I+m*l^2)/p; 0; m*l/p];

C = [1 0 0 0]; % apenas o angulo do pendulo é medido

D = 0;

sys_ss_c = ss(A,B,C,D)

% Modelo de Espaço de Estado Discreto
Ts = 0.01; % s
sys_ss_d = c2d(sys_ss_c,Ts)
[Ad,Bd,Cd,Dd] = ssdata(sys_ss_d);

% Modelo Aumentando de Espaço de Estado Discreto
% Adicionar o integrador

% Op. 01: Integrar a saída 
Aa = [1,     Cd*Ad;
       zeros(size(Ad,1),1), Ad];  
Ba = [Cd*Bd;
       Bd];
Ca = [1, zeros(1, size(Ad,1))];     % só mede a variavel integrada
Da = 0;

% Op. 02: Integrar o erro
% Aa = [Ad,          zeros(size(Ad,1),1);
%      -Cd*Ad,       1       ];  % erro = r - y = r - Cd*x
% Ba = [Bd;
%      -Cd*Bd];
% Ca = [Cd, 0];

sys_ss_d_a = ss(Aa,Ba,Ca,Da,Ts)
%% LQR

% Análise de estabilidade relativa via margens de ganho e de fase

% Teste de seguimentos de referência e saída medida

%% Filtro de Kalman

% Análise de estabilidade relativa via margens de ganho e de fase

% Teste de seguimentos de referência e saída medida



%% Compensador dinâmico completo LQG

% Análise de estabilidade relativa via margens de ganho e de fase



