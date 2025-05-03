%% Atividade A1: Estudo de caso para o controle de velocidade longitudinal de um quadricóptero.
clear all; close all; clc;

%% Setup
% Carregando os coeficientes do modelo de 2ª ordem
load('sys_mod01.mat');
Ts = model.Ts;
Az = model.Az; a1 = Az(2); a2 = Az(3);
Bz = model.Bz; b0 = Bz(2); b1 = Bz(3);
%Gz = tf(Bz,Az,Ts);

% Configurações
tfinal = 60; % em segundos
N = round( tfinal/Ts ); % numero total de amostras

%% Controlador PID
% Utilizando o modelo de 2ª ordem para sintonia 
kp = -1.2; ki = -0.5; kd =-0.01;

% Aproximação Backward Difference
s0 = kp +ki*Ts +kd/Ts;
s1 = -kp -2*kd/Ts;
s2 = kd/Ts;

% Sinal de referencia (degrau unitario)
r = zeros(N,1); r(round(N/4):end)=1; 

% Condições iniciais e inicializacao
ym = zeros(N,1); um = zeros(N,1);
for k = 1:2
    r(k) = 0; ym(k) = 0; um(k) = 0;  
end

% Teste
for k = 3:N
    % Planta
    ym(k) = -a1*ym(k-1) -a2*ym(k-2) +b0*um(k-1) +b1*um(k-2); 

    %% Calcular o sinal de controle 
    um(k) = um(k-1) + s0*(r(k)-ym(k)) + s1*(r(k-1)-ym(k-1)) + s2*(r(k-2)-ym(k-2));
    
    % Limitando o sinal de controle 
    if um(k) <= -1
        um(k) = -1;
    elseif um(k) >= 1
        um(k) = 1;
    end

end

% Plot
t = (0:N-1)' * Ts;

subplot(2,1,1)
plot(t, r, 'k--', 'LineWidth', 1.2); hold on;
plot(t, ym, 'b', 'LineWidth', 1.5);
legend('Referência r(k)', 'Saída ym(k)', 'Location', 'northeast');
xlabel('Tempo [s]');
ylabel('Velocidade [m/s]');
title('Resposta do Sistema ao Degrau Unitário');
grid on;

subplot(2,1,2)
plot(t, um, 'r', 'LineWidth', 1.5);
xlabel('Tempo [s]');
ylabel('Sinal de Controle u(k)');
title('Ação de Controle PID');
grid on;

%% Analise das restricoes
% Sobressinal máximo
overshoot = ((max(ym) - r(end)) /  r(end)) * 100;
fprintf('Overshoot: %.2f %%\n', overshoot);

% Tempo de assentamento
tol = 0.05 * r(end); % tolerancia de 5%

idx_start = find(r > 0, 1, 'first');  % indice do inicio do degrau
idx_settle = find(abs(ym - r(end)) > tol, 1, 'last'); % indice onde y se aproxima de r

t_settle = t(idx_settle) - t(idx_start);
fprintf('Settling time: %.2f s\n', t_settle);

%% Salvando o controlador
control.kp = kp;
control.ki = ki;
control.kd = kd;
%control.Cz = Cz;
control.Ts = Ts;
save 'control_pid.mat' control;