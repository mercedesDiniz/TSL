%% Atividade 4: Controle LQG/LTR para a velocidade lateral e longitudinal de um quadrirotor.
clear all; close all; clc;

%% Planta do sistema do tipo quadrirotor com saturação nos atuadores 
% Ref.: https://lacos.ufpa.br/plantas/ardrone_sim/ardrone_sim.html

path_my = 'D:\PPGEE-UFPA\2025.1-TSL-Antonio\'; 
path_lasse = 'C:\Users\Mercedes Diniz\Documents\PPGEE-UFPA\TSL-2025.1[Silveira]\';
path_tsl = 'atividades\a4\modelo\quadrotor_model.mat';

load([path_my, path_tsl]);

% x1: phi (rad): ângulo de rolagem
% x2: theta (rad): ângulo de inclinação
% x3: u_spd (m/s): velocidade longitudinal
% x4: v_spd (m/s): velocidade lateral

% u1: thrust lateral: saturações [-1;1]
% u2: thrust longitudinal: saturações [-1;1]

%% Modelo em Espaço de Estado em Tempo Discreto
Ts = quadrotor_model.Ts

A = quadrotor_model.A(1:4,1:4);
B = quadrotor_model.B(1:4,1:2);
C = [0 0 0 1;  % y1 = x4: v_spd (m/s): velocidade lateral
     0 0 1 0]; % y2 = x3: u_spd (m/s): velocidade longitudinal
D = zeros(2,2);

sys = ss(A,B,C,D,Ts);

n = length(A);  % ordem do sistema
ny = 2;         % número de saídas
nu = 2;         % número de entradas

%% Modelo Aumentando em Espaço de Estado Discreto (adicionar o integrador e pré-compensador DC)
% P = inv( dcgain(sys) )
P = inv( C*inv(eye(n,n) -A)*B ); % Pré-compensador

Aa = [  eye(ny,ny)     C*A;
       zeros(4,ny)    A  ];

n = length(Aa);  % ordem do sistema

Ba = [C*B; B]*P;

Ca = [eye(ny,ny) zeros(ny,4)];

Da = D;

sysa = ss(Aa,Ba,Ca,Da,Ts);

figure; sigma(sysa);

%% Verificação da Controlabilidade e Observabilidade
Co = ctrb(Aa, Ba);  % Matriz de Controlabilidade
Ob = obsv(Aa, Ca);  % Matriz de Observabilidade

disp('Verificação da Controlabilidade:');
disp(['> rank(Co) = ' num2str(rank(Co))]);
if rank(Co) == n
    disp(['== ' num2str(n) ' - O sistema é controlável' ]);
else
    disp(['!= ' num2str(n) ' - O sistema não é controlável' ]);
end

disp('Verificação da Observabilidade:');
disp(['> rank(Ob) = ' num2str(rank(Ob))]);
if rank(Ob) == n
    disp(['== ' num2str(n) ' - O sistema é observável' ]);
else
    disp(['!= ' num2str(n) ' - O sistema não é observável' ]);
end

% %% Controle via realimentação de estado estimado - LQR
% disp('Controlador LQR');
% Q_lqr = diag([ 1   1    1     1      1 ]);  % peso das variaveis 
% R_lqr = 1;                                  % peso do esforço de controle
% 
% %% Estimador de estados - Filtro de Kalman
% disp('Filtro de Kalman');
% Q_kf = diag([ 1   1    1     1      1 ]);  % peso das variaveis 
% R_kf = 1;                                  % ruído da medição
% 
% 
% %% Simulação do Compensador dinâmico completo - LQG
% t_f = 20;           % tempo de simulação
% N = round(t_f/Ts);  % número de amostras
% t = 0:Ts:N*Ts-Ts;   % vetor de tempo discreto



