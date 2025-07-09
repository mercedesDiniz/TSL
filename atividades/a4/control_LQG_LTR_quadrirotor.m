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

Ba = [C*B; B]*P;

Ca = [eye(ny,ny) zeros(ny,4)];

Da = D;

sysa = ss(Aa,Ba,Ca,Da,Ts);

figure; sigma(sysa);
title('Principais ganhos (ou valores singulares) do sistema MIMO');

%% Verificação da Controlabilidade e Observabilidade da realização
na = length(Aa);    % ordem do sistema aumentado
Co = ctrb(Aa, Ba);  % Matriz de Controlabilidade
Ob = obsv(Aa, Ca);  % Matriz de Observabilidade

disp('Verificação da Controlabilidade:');
disp(['> rank(Co) = ' num2str(rank(Co))]);
if rank(Co) == na
    disp(['== ' num2str(na) ' - O sistema é controlável' ]);
else
    disp(['!= ' num2str(na) ' - O sistema não é controlável' ]);
end

disp('Verificação da Observabilidade:');
disp(['> rank(Ob) = ' num2str(rank(Ob))]);
if rank(Ob) == na
    disp(['== ' num2str(na) ' - O sistema é observável' ]);
else
    disp(['!= ' num2str(n) ' - O sistema não é observável' ]);
end

%% Projeto do filtro de Kalman
    % x1: phi (rad): ângulo de rolagem
    % x2: theta (rad): ângulo de inclinação
    % x3: u_spd (m/s): velocidade longitudinal
    % x4: v_spd (m/s): velocidade lateral 
    % y1 = x4: v_spd (m/s): velocidade lateral
    % y2 = x3: u_spd (m/s): velocidade longitudinal

%            y1    y2   dx1 dx2 dx3 dx4
Q_kf = diag([ 1     1    1   1   1   1 ]);

%            y1     y2
R_kf = diag([ 1e5   1e5 ]);

L = ( dlqr(Aa',Ca',Q_kf, R_kf) )'

% Analise do filtro de Kalman em malha fechada
sys_kf = ss( Aa-L*Ca, L, Ca, Da, Ts );
    % Avaliação das propriedades de rastreamento no domínio do tempo
    figure; step(sys_kf);
    title('Propriedades de rastreamento do filtro de Kalman');
    
    % Avaliação do sistema no domínio da frequência
    Tsen_kf = sys_kf;
    Ssen_kf = eye(2,2) -Tsen_kf;
    figure; sigma(Tsen_kf); hold; sigma(Ssen_kf); grid;
    title('Funções de sensibilidade do filtro de Kalman');
    legend('|T(e^{j\omegaT_s})|','|S(e^{j\omegaT_s})|');
    
    % Picos de sensibilidade
    mt_kf = max( max( sigma(Tsen_kf) ) );
    ms_kf = max( max( sigma(Ssen_kf) ) );
    
    % Margens de ganho e fase
    GmdB_kf = min( 20*log10(ms_kf/(ms_kf-1)), 20*log10(1+(1/mt_kf)) );
    Pmdeg_kf = (180/pi)*min( (2*asin(1/(2*ms_kf)) ), (2*asin(1/(2*mt_kf)) ) );

    disp('Margens de ganho e fase do filtro de Kalman:');
    disp('GmdB = '); disp(GmdB_kf);
    disp('Pmdeg = '); disp(Pmdeg_kf);

