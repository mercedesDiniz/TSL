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
Ts = quadrotor_model.Ts;

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
Q_kf = diag([ 1     1    1   1   1e-3   1e-3 ]);

%            y1     y2
R_kf = diag([ 1e4   1e4 ]);

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
    mt = max( max( sigma(Tsen_kf) ) );
    ms = max( max( sigma(Ssen_kf) ) );
    
    % Margens de ganho e fase
    GmdB_kf = min( 20*log10(ms/(ms-1)), 20*log10(1+(1/mt)) );
    Pmdeg_kf = (180/pi)*min( (2*asin(1/(2*ms)) ), (2*asin(1/(2*mt)) ) );

    disp('Margens de ganho e fase do filtro de Kalman:');
    disp('GmdB = '); disp(GmdB_kf);
    disp('Pmdeg = '); disp(Pmdeg_kf);

%% Projeto do LQR
    % x1: phi (rad): ângulo de rolagem
    % x2: theta (rad): ângulo de inclinação
    % x3: u_spd (m/s): velocidade longitudinal
    % x4: v_spd (m/s): velocidade lateral 
    % y1 = x4: v_spd (m/s): velocidade lateral
    % y2 = x3: u_spd (m/s): velocidade longitudinal

 %            y1    y2   dx1 dx2 dx3 dx4
Qlq = diag([ 1     1    1   1   1   1 ]);

%            y1  y2
Rlq = diag([ 1   1 ]);


K = dlqr(Aa,Ba,Qlq,Rlq)

% Analise do filtro do LQR em malha fechada
sys_lqr = ss(Aa-Ba*K, Ba*K(1:2,1:2), Ca, Da, Ts);
    % Avaliação das propriedades de rastreamento no domínio do tempo
    figure; step(sys_lqr);
    title('Propriedades de rastreamento do LQR');
    
    % Avaliação do sistema no domínio da frequência
    Tsen_lqr = sys_lqr;
    Ssen_lqr = eye(2,2) - Tsen_lqr;
    figure; sigma(Tsen_lqr); hold; sigma(Ssen_lqr); grid;
    title('Funções de sensibilidade do LQR');
    legend('|T(e^{j\omegaT_s})|','|S(e^{j\omegaT_s})|');
    
    % Picos de sensibilidade
    mt = max( max( sigma(Tsen_lqr) ) );
    ms = max( max( sigma(Ssen_lqr) ) );
    
    % Margens de ganho e fase
    GmdB_lqr = min( 20*log10(ms/(ms-1)), 20*log10(1+(1/mt)) );
    Pmdeg_lqr = (180/pi)*min( (2*asin(1/(2*ms)) ), (2*asin(1/(2*mt)) ) );

    disp('Margens de ganho e fase do LQR:');
    disp('GmdB = '); disp(GmdB_lqr);
    disp('Pmdeg = '); disp(Pmdeg_lqr);

%% Simulação do LQG
t_f = 20;           % tempos de simulação
N = round(t_f/Ts);  % número de amostras
t = 0:Ts:N*Ts-Ts;   % vetor de tempo discreto

    % Sinal de referencia
    r1 = zeros(1,N); r1(round(N/3):end) = 1;   % v_spd (x4)
    r2 = zeros(1,N); r2(round(N/3):end) = 1;   % u_spd (x3)
    r = [r1; r2];                       

    % Disturbios na entrada e saida
    w = 0 * randn(n, N);  
    v = 0 * randn(ny, N);

    % Condições iniciais do modelo nominal
    x = zeros(n, N);
    y = zeros(ny, N);
    u = zeros(nu, N);
    du = zeros(nu, N); 

    % Condições iniciais do modelo aumentado
    xa = zeros(na, N); xa(:,1) = zeros(na,1);
    ya = zeros(ny, N); 

    % Condições iniciais das variaveis estimadas
  
for k = 2:N
    % Modelo em espaço de estados
    x(:,k)     = A * x(:,k-1) + B * du(:,k-1) + w(:,k-1);
    y(:,k-1)   = C * x(:,k-1) + v(:,k-1); 

    % Modelo em espaço de estados aumentado com o filtro de Kalman
    ya(:,k-1) = Ca * xa(:,k-1);
    xa(:,k) = Aa * xa(:,k-1) + Ba * du(:,k-1) + L * (r(:,k-1) - ya(:,k-1));

    % Lei de controle de realimentação de estados
     du(:,k) = K(:,1:ny) * (r(:,k) - xa(1:ny, k));

    % Passando du(k) pelo integrador discreto 
    u(:,k) = u(:,k-1) + du(:,k);

    % Saturação dos atuadores
    u(:,k) = min(max(u(:,k), -1), 1);
end

%% Plots

figure;
subplot(3,1,1);
plot(t, r1, 'k--', t, y(1,:), 'b', t, ya(1,:), ':r');
ylabel('v_{spd} (m/s)');
legend('Ref', 'Saída', 'Estimada');
grid on;

subplot(3,1,2);
plot(t, r2, 'k--', t, y(2,:), 'b', t, ya(2,:), ':r');
ylabel('u_{spd} (m/s)');
legend('Ref', 'Saída', 'Estimada');
grid on;

subplot(3,1,3);
plot(t, u(1,:), 'r', t, u(2,:), 'b');
xlabel('Tempo (s)');
ylabel('Controle');
legend('u_1 (lat)', 'u_2 (long)');
grid on;
sgtitle('Simulação LQG MIMO');


