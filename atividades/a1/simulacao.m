%% Atividade A1: Estudo de caso para o controle de velocidade longitudinal de um quadricóptero.
clear all; close all; clc;

%% Setup
% Carregando os coeficientes do modelo de 3ª ordem
load('sys_mod02.mat');
Ts = model.Ts;
Az = model.Az; a1 = Az(2); a2 = Az(3); a3 = Az(4);
Bz = model.Bz; b0 = Bz(2); b1 = Bz(3); b2 = Bz(4); 

% Carregando os valores de kp, ki e kd
load('control_pid.mat');
kp = control.kp;  ki = control.ki; kd = control.kd;

% Aproximação Backward Difference
s0 = kp +ki*Ts +kd/Ts; s1 = -kp -2*kd/Ts; s2 = kd/Ts;

% Configurações
tfinal = 60; % em segundos
N = round( tfinal/Ts ); % numero total de amostras

% Sinal de referencia (degrau unitario)
r = zeros(N,1); r(round(N/4):end)=1; 

% Condições iniciais e inicializacao
ym1 = zeros(N,1); ym2 = zeros(N,1);
um1 = zeros(N,1); um2 = zeros(N,1);
for k = 1:3
    r(k) = 0; ym1(k) = 0; um1(k) = 0;  
end

%% Perturbação de carga: 30% da referência
pert = zeros(N,1);
%pert(round(30/Ts):end) = 0.3 * r(end);

%% Simulacao
for k = 4:N
    % Planta
    ym1(k) = -a1*ym1(k-1) -a2*ym1(k-2) -a3*ym1(k-3) +b0*um1(k-1) +b1*um1(k-2) +b2*um1(k-3) + pert(k);  
    ym2(k) = -a1*ym2(k-1) -a2*ym2(k-2) -a3*ym2(k-3) +b0*um2(k-1) +b1*um2(k-2) +b2*um2(k-3) + pert(k);

    % Calcular o sinal de controle do PID
    um1(k) = um1(k-1) + s0*(r(k)-ym1(k)) + s1*(r(k-1)-ym1(k-1)) + s2*(r(k-2)-ym1(k-2));

    % Calculo do sinal de controle do PI-D
    um2(k) = um2(k-1) +(kp + ki*Ts)*r(k) -kp*r(k-1) - s0*ym2(k) - s1*ym2(k-1) - s2*ym2(k-2);
    
    % Limitando o sinal de controle 
    if um1(k) <= -1
        um1(k) = -1;
    elseif um1(k) >= 1
        um1(k) = 1;
    end

    if um2(k) <= -1
        um2(k) = -1;
    elseif um2(k) >= 1
        um2(k) = 1;
    end
end

%% Analise do desempenho
% Erro
erro1 = r - ym1; erro2 = r - ym2;

% Integral of Squared Error 
ISE1 = sum(erro1.^2);    ISE2 = sum(erro2.^2);

% Normalized Mean Squared Error
NMSE1 = sum(erro1.^2) / sum(r.^2);  NMSE2 = sum(erro2.^2) / sum(r.^2);

% Energia e potência do sinal de controle
e_u1 = sum(um1.^2); e_u2 = sum(um2.^2);
pot_u1 = e_u1 / N;  pot_u2 = e_u2 / N;

fprintf('\n=== Avaliação dos Índices ===\n');
fprintf('Controlador PID:    ISE = %.2f, NMSE = %.2f, Energia = %.2f, Potência = %.2f\n', ISE1, NMSE1, e_u1, pot_u1);
fprintf('Controlador PI-D:   ISE = %.2f, NMSE = %.2f, Energia = %.2f, Potência = %.2f\n', ISE2, NMSE2, e_u2, pot_u2);

%% Plots
figure;
subplot(2,1,1);
plot((0:N-1)*Ts, r, '--k', (0:N-1)*Ts, ym1, 'b', (0:N-1)*Ts, ym2, 'r');
legend('Referência', 'y - PID', 'y - PI-D');
ylabel('Amplitide (V)'); title('Resposta dos Sistemas');

subplot(2,1,2);
plot((0:N-1)*Ts, um1, 'b', (0:N-1)*Ts, um2, 'r');
legend('PID', 'PI-D');
xlabel('Tempo (s)'); ylabel('u[k]'); title('Sinais de Controle');

