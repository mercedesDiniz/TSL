%% Atividade 02: Controle do pêndulo invertido sobre carro usando LQG.
clear all; close all; clc;

%% Planta do pêndulo invertido 
% Ref.: https://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling 
% A variável sensorialmente medida é a posição angular do pêndulo (em rad).
% Sinal de controle (u) é a Força aplicada sobre o carro (F) em Newtons.
% As variaveis de estado são:
    % x1 = x:       posição do carro em metros
    % x2 = x_dot:   velocidade do carro em m/s
    % x3 = phi:     angulo do pendulo (em rad)
    % x4 = phi_dot  velocidade angular do pendulo (rad/s)

%% Modelo em Espaço de Estados em Tempo Continuo
M = 0.5; m = 0.2; b = 0.1; I = 0.006; g = 9.8; l = 0.3; 
p = I*(M+m)+M*m*l^2; 

% Matrizes de estado
A = [0      1              0           0;
     0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
     0      0              0           1;
     0 -(m*l*b)/p       m*g*l*(M+m)/p  0];

B = [   0; 
    (I+m*l^2)/p;
        0; 
      m*l/p];

C = [0 0 1 0]; ny = size(C,1); % sersor da posição angular do pêndulo
% C = [0 0 1 0; 0 0 1 0]; ny = size(C,1); % + sensor de posição do carro

D = zeros(ny,1);
    
sys_ss_c = ss(A,B,C,D, ...
    'statename', {'x' 'x_dot' 'phi' 'phi_dot'}, ...
    'inputname',{'u'}, ...
    'outputname', {'phi'})

% Analise dos poles e autovalores do sistema
disp('Autovalores contínuos do sistema:');
disp(eig(A)); % o polo 5.5651 é instavel (está no semiplano direito do plano s)

% Plot dos polos no plano S
figure; rlocus(sys_ss_c);
title('Disposição dos polos do modelo contínuo'); grid on; hold on;
theta = linspace(0, 2*pi, 100);
plot(cos(theta), sin(theta), '--', 'Color', [0.5 0.5 0.5]);
axis equal;
legend('Lugar da raiz', 'Círculo Unitário');

% Analise das maximas frequencias do sistema
disp('Verificando a frequência máxima do sistema:');
damp(sys_ss_c) 
f_max = max(damp(sys_ss_c));
disp(['> f_max = ' num2str(f_max) ' rad/s']);

%% Modelo em Espaço de Estado em Tempo Discreto 
v_factor = 200;                  % fator de velocidade da amostragem  
fs = (f_max*v_factor)/(2*pi);    % frequencia de amostragem (Hz) 
Ts = 1/fs;                       % periodo de amostragem (s) 
Ts = fix(v_factor*Ts)/v_factor; 

disp(['Período de Amostragem Ts = ' num2str(Ts) 's']);

% Discretizando o sistema
sys_ss_d = c2d(sys_ss_c, Ts)
[Ad, Bd, Cd, Dd] = ssdata(sys_ss_d);

% Compara o modelo contínuo com o discreto 
% Observa-se que os modelos divergem nas altas frequências.
% O sistema é instavel, pq na região de inversão de fase há o aumento do ganho.
figure; bode(sys_ss_c, sys_ss_d); title('Digrama de bode do modelo cont. e disc.');
legend('Contínuo', 'Discreto', 'Location','southwest');

%% Modelo Aumentando em Espaço de Estado Discreto (adicionar o integrador)
Aa = [eye(ny, ny)   Cd*Ad;
      zeros(4, ny)  Ad];

Ba = [Cd*Bd;
       Bd];

Ca = [eye(ny, ny) zeros(ny, 4)];  

Da = Dd;

sys_ss_d_a = ss(Aa, Ba, Ca, Da, Ts,  ...
    'statename', {'x' 'x_dot' 'phi' 'phi_dot' 'int'}, ...
    'inputname',{'u'}, ...
    'outputname', {'phi'})

figure; bode(sys_ss_d_a); grid; title('Digrama de bode do modelo aumentado');

% Analise dos poles e autovalores do sistema
disp('Autovalores discretos do sistema aumentado:');
disp(eig(Aa)); % o polo 1.0282 é instavel (está fora do circulo únitario no plano z)
               % aparentemente a variavel na posição pendulo já é uma
               % variavel integradora (a mesma se acumula ao longo do tempo)

%% Testes de Controlabilidade e Observabilidade
n = length(Aa);     % Número de variaveis de estado (o integrador é considerado)
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

%% Controlador LQR
disp('Controlador LQR');
             % y1  x  x_dot  phi  phi_dot  
Q_lqr = diag([ 1   1    1     1      1 ]);  % peso das variaveis 
R_lqr = 1;                                  % peso do esforço de controle

disp('Verificação da Detectabilidade (p/ o LQR):');
Ob_dec_lqr = obsv(Aa', sqrt(Q_lqr)');
disp(['> rank(obsv(Aa^(T), sqrt(Q_lqr)^(T))) = ' num2str(rank(Ob_dec_lqr))]);

if rank(Ob_dec_lqr) == n
    disp(['== ' num2str(n) ' - O sistema é detectável' ]);
else
    disp(['!= ' num2str(n) ' - O sistema não é detectável' ]);
end

% Calculo da matriz P pela Equação a Diferenças de Riccati
P = 1e6*eye(n, n); p_samples = 1000; 
for i = 1:p_samples
    P = Aa'*P*Aa -Aa'*P*Ba*inv(Ba'*P*Ba+R_lqr)*Ba'*P*Aa+Q_lqr;
    traceP(i)=trace(P); 
end
figure; plot(traceP); ylabel('Traço de P - LQR'); xlabel('Amostras');

% Calculo do Ganho
K = ( Aa'*P*Ba*inv(Ba'*P*Ba+R_lqr) )';
disp('LQG K = '); disp(K);

% Analise dos poles e autovalores do sistema
disp('Autovalores do controlador LQR:');
disp(eig(Aa-Ba*K)); % precisam ser estaveis (está dentro do circulo únitario no plano z)

% TODO: Analise em MF pelas função de sensibilidade
% Tsen = ss(Aa-Ba*K,Ba*K(4),Ca,Da,Ts);
% Ssen = eye(ny,ny)-Tsen; 
% 
% mt = max( sigma(Tsen) );
% ms = max( sigma(Ssen) );
% 
% GmdB = min( 20*log10(ms/(ms-1)), 20*log10(1+(1/mt)) );
% Pmdeg = (180/pi)*min( (2*asin(1/(2*ms)) ), (2*asin(1/(2*mt)) ) );
% 
% disp('Margens de ganho e fase obtidas por análise de malha fechada:');
% disp('GmdB = '); disp(GmdB);
% disp('Pmdeg = '); disp(Pmdeg);
% 
% figure; sigma(Tsen); hold; sigma(Ssen); grid;
% legend('|Tsen|','|Ssen|'); title('LQR Sensibilidades');

%% Filtro de Kalman (Observador)
disp('Filtro de Kalman');
             % y1  x  x_dot  phi  phi_dot  
Q_kf = diag([ 1   1    1     1      1 ]);  % peso das variaveis 
R_kf = 1;                                  % ruído da medição

disp('Verificação da Detectabilidade (p/ o Filtro de Kalman):');
Ob_dec_kf = obsv(Aa', sqrt(Q_kf)');
disp(['> rank(obsv(Aa^(T), sqrt(Q_kf)^(T))) = ' num2str(rank(Ob_dec_kf))]);

if rank(Ob_dec_kf) == n
    disp(['== ' num2str(n) ' - O sistema é detectável' ]);
else
    disp(['!= ' num2str(n) ' - O sistema não é detectável' ]);
end

% Calculo da matriz P pela Equação a Diferenças de Riccati
Po = 1e6*eye(n, n); p_samples = 1000; 
for i = 1:p_samples
    Po = Aa*Po*Aa' -Aa*Po*Ca'*inv(Ca*Po*Ca'+R_kf)*Ca*Po*Aa'+Q_kf;
    tracePo(i)=trace(Po); 
end
figure; plot(tracePo); ylabel('Traço de P - Filtros de Kalman'); xlabel('Amostras');

% Calculo do Ganho
L = ( Aa*Po*Ca'*inv(Ca*Po*Ca'+R_kf) );
disp('KF L = '); disp(L);

% Analise dos poles e autovalores do sistema
disp('Autovalores do filtros de Kalman:');
disp(eig(Aa-L*Ca)); % precisam ser estaveis (está dentro do circulo únitario no plano z)

% TODO: Analise em MF pelas função de sensibilidade

%% Simulação do controlador LQG
t_f = 20;           % temos de simulação
N = round(t_f/Ts);  % número de amostras
t = 0:Ts:N*Ts-Ts;   % vetor de tempo discreto
    
    % Sinal de referencia
    r1(1:5) = 0; r1(6:N) = 0;

    % Disturbios na entrada e saida
    w1 = 0*wgn(1, N, 1e-3, 'linear'); % ruido
    w2(1:round(N/2)) = 0; w2(round(N/2)+1:N) = 0*0.02;

    v1 = 0*wgn(1, N, 1e-3, 'linear'); % ruido
    v2(1:round(N/2)) = 0; v2(round(N/2)+1:N) =+0.2*0;

    % Condições iniciais do modelo nominal
                % x  x_dot    phi    phi_dot
    x(:,:, 1) = [0;  0;    0.0349;     0];   % entrada
    y(:,:, 1) = C*x(:,:, 1);                 % saida
    u(1) = 0;                                % sinal de controle

    % Condições iniciais do modelo aumentado
                % y1  x  x_dot  phi  phi_dot
    xa(:,:, 1) = [0;  0;    0;     0;     0];   % entrada
    ya(:,:, 1) = Ca*xa(:,:, 1);                 % saida
    du(1) = 0;                                  % sinal de controle

    % Condições iniciais das variaveis estimadas
    x1e(1)= 0; x2e(1)= 0; x3e(1)= 0; x4e(1)= 0;

for k = 2:N
    % Modelo em espaço de estados
    x(:,:, k) = A*x(:,:, k-1) + B*u(k-1) + [1; 0; 0; 0]*w2(k-1); 
    y(:,:, k) = C*x(:,:, k) + v1(k) + v2(k);

    % Modelo em espaço de estados aumentado com o filtro de Kalman
    xa(:,:, k) = Aa*xa(:,:, k-1) + Ba*du(k-1) + L*(y(:,:, k-1) - ya(:,:, k-1));
    ya(:,:, k) = Ca*xa(:,:, k);

    % Estimação das variaveis
    x1e(k) = x1e(k-1) + xa(2, k);
    x2e(k) = x2e(k-1) + xa(3, k);
    x3e(k) = x3e(k-1) + xa(4, k);
    x4e(k) = x4e(k-1) + xa(5, k);

    % Lei de controle de realimentação de estados
    du(k) = K*([r1(k); 0; 0; 0; 0] - xa(:,:, k));

    % Passando du(k) pelo integrador discreto   
    u(k) = u(k-1) + du(k);
end

% Plots
figure;
subplot(211)
    plot(t, r1, 'k'); hold;
    plot(t, y(1,:),'b'); plot(t, ya(1,:), ':r');
    legend('r1(t)', 'y1(t)', 'y1a(t): estimado');
    ylabel('Angulo (rad)');
subplot(212)
    plot(t, u, 'b'); xlabel('Tempo (s)'); ylabel('Força (N)');

figure;
subplot(411)
    plot(t,xa(1,:),'r');
    ylabel('xa_1(t)');
    title('Variáveis ​​de estado do estimador aumentado');

subplot(412)
    plot(t,xa(2,:),'r');
    ylabel('xa_2(t)');

subplot(413)
    plot(t,xa(3,:),'r');
    ylabel('xa_3(t)');

subplot(414)
    plot(t,xa(4,:),'r');
    ylabel('xa_4(t)');
