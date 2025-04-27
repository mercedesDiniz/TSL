%% Atividade A1: Estudo de caso para o controle de velocidade longitudinal de um quadricóptero.
clear all; close all; clc;

%% Setup
% Carregando o datalog (dados reais)
data = load('datalog_lon_vel.txt');
t = data(:,1);          % Coluna 1: tempo (s)
u = data(:,2);          % Coluna 2: comando [-1,1]
y = data(:,3);          % Coluna 3: velocidade [m/s]

% Carregando os controladores
load('control_pid.mat');
Cz1 = control.Cz; 

load('control_pi-d.mat');
Cz2 = control.Cz;

% Carregando o modelo de 3ª ordem da planta
load('sys_mod02.mat');
Ts = model.Ts;
Az = model.Az; a1 = Az(2); a2 = Az(3); a3 = Az(4);
Bz = model.Bz; b0 = Bz(2); b1 = Bz(3); b2 = Bz(4);

Gz = tf(Bz,Az,Ts);

% Configurações da simulação
tfinal = 20; % em segundos
N = round( tfinal/Ts ); % numero total de amostras

%% Analise do sistema com o controlador PID
disp('SYSTEM WITH PID CONTROLLER')
% Analise no dominio da frequencia
figure; margin(Cz1*Gz); 

% Funcao de transferencia de malha fechada
Gmfz1 = feedback(Cz1*Gz,1,-1); % C(z)*G(z) / 1-C(z)*G(z) 
disp('Closed-loop poles: '); pole(Gmfz1)
disp('Closed-loop zeros: '); zero(Gmfz1)

% Analise de estabilidade relativa em malha fechada
% Calculo das margens de ganho e de fase aproximadas
Tsen = Gmfz1;           % função de co-sensibilidade
Ssen = 1 -Tsen;         % função de sensibilidade
mt = max( sigma(Tsen) ); ms = max( sigma(Ssen) );
GmdB = min( 20*log10(ms/(ms-1)), 20*log10(1+(1/mt)) );
Pmdeg = (180/pi)*min( (2*asin(1/(2*ms)) ), (2*asin(1/(2*mt)) ) );

disp('Gain and Phase margins obtained by closed-loop analysis:');
disp('GmdB = '); disp(GmdB);
disp('Pmdeg = '); disp(Pmdeg);

% Plot da func. de sensibilidade e co-sensibilidade no dominio da freq
%figure; sigma(Tsen); hold; sigma(Ssen); grid;
%legend('|Tsen|','|Ssen|');

%% Simulacao do sistema em malha fechada com o controlador PID
% Gerando o sinal de referencia

%% Analise do sistema com o controlador PI-D

%% Simulacao do sistema em malha fechada com o controlador PI-D

