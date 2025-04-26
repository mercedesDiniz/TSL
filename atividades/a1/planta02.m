%% Atividade A1: Estudo de caso para o controle de velocidade longitudinal de um quadricóptero.
%clear all; close all; clc;

%% Setup
% Carregando os dados
data = load('datalog_lon_vel.txt');
t = data(:,1);          % Coluna 1: tempo (s)
u = data(:,2);          % Coluna 2: comando [-1,1]
y = data(:,3);          % Coluna 3: velocidade [m/s]

t_sim = 60;             % tempo de simulação
Ts = 0.065;             % periodo de amostragem
N = round( t_sim/Ts ); %length(t)    % numero total de amostras

% Plot do Datalog
figure; 
% Velocidade longitudinal
subplot(2,1,1); plot(t, y, 'b'); grid on;
xlabel('Time [s]'); ylabel('Longitudinal Velocity [m/s]'); 
title('DATALOG', 'FontWeight', 'bold');
% Comando
subplot(2,1,2); plot(t, u, 'b'); grid on;
xlabel('Time [s]'); ylabel('Longitudinal Thrust')

%% Modelagem da Planta - 3ª Ordem
% G(z) = Y(Z)/U(Z) = b0*z^(-1) + b1*z^(-2) + b2*z^(-3) / 1 + a1*z^(-1) + a2*z^(-2) + a3*z^(-3)
% y[n] = -a1*ym(n-1) -a2*y(n-2) -a3*y(n-3) + b0*u(n-1) + b1*u(n-2) + b2*u(n-3);

% Estimação dos parametros da eq. das diferenças

% Fatiando o datalog
t_i = 10; N_i = round(t_i/Ts);
t_f = 20; N_f = round(t_f/Ts);
u_i = u(N_i:N_f); y_i = y(N_i:N_f); N_id = length(y_i);

% Estimador dos Minimos Quadrados 
PHI=[];
for n=4:N_id
    PHI = [PHI ; [-y_i(n-1) -y_i(n-2) -y_i(n-3) u_i(n-1) u_i(n-2) u_i(n-3)] ];
end

theta = inv(PHI' * PHI) * PHI' * y_i(4:N_id);
a1 = theta(1); a2 = theta(2); a3 = theta(3); b0 = theta(4); b1 = theta(5); b2 = theta(6);

%% Simulação da saída estimada pelo modelo
ym = zeros(N,1);          % inicialização da saida estimada
ym(1:3) = y(1:3);         % condições iniciais

for n = 4:N
    ym(n) = -a1*ym(n-1) -a2*ym(n-2) -a3*ym(n-3) +b0*u(n-1) +b1*u(n-2) +b2*u(n-3);
end

%% Plot do modelo identificado
figure;
%plot(t, u, 'r', 'LineWidth', 1); hold on;
plot(t, y, 'b', 'LineWidth', 1);  hold on;
plot(t(1:length(t)-1), ym, 'k', 'LineWidth', 1);

legend('y(t)', 'ym(t)', 'Location', 'northeast');
xlabel('Time (s)'); ylabel('Amplitude');
title('Comparação entre o Sistema Real e o Modelo de 3ª Ordem');
grid on;



