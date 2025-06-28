%% Atividade 03: Análise da estabilidade relativa do controle digital do pêndulo invertido sobre carro (A2).
clear all; close all; clc;

%% Planta do pêndulo invertido 
% Ref.: https://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling 

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

D = zeros(ny,1);
    
sys_ss_c = ss(A,B,C,D, ...
    'statename', {'x' 'x_dot' 'phi' 'phi_dot'}, ...
    'inputname',{'u'}, ...
    'outputname', {'phi'})

f_max = max(damp(sys_ss_c));


%% Modelo em Espaço de Estado em Tempo Discreto 
v_factor = 200;                  % fator de velocidade da amostragem  
fs = (f_max*v_factor)/(2*pi);    % frequencia de amostragem (Hz) 
Ts = 1/fs;                       % periodo de amostragem (s) 
Ts = fix(v_factor*Ts)/v_factor; 

% Discretizando o sistema
sys_ss_d = c2d(sys_ss_c, Ts)
[Ad, Bd, Cd, Dd] = ssdata(sys_ss_d);

%% Modelo Aumentando em Espaço de Estado Discreto (adicionar o integrador)
Aa = [eye(ny, ny)   Cd*Ad;
      zeros(4, ny)  Ad];  n = length(Aa);     % Número de variaveis de estado (o integrador é considerado)

Ba = [Cd*Bd;
       Bd];

Ca = [eye(ny, ny) zeros(ny, 4)];  

Da = Dd;

sys_ss_d_a = ss(Aa, Ba, Ca, Da, Ts,  ...
    'statename', {'x' 'x_dot' 'phi' 'phi_dot' 'int'}, ...
    'inputname',{'u'}, ...
    'outputname', {'phi'})

%% Controle via realimentação de estado estimado - LQR
disp('Controlador LQR');
%Teste 0
             % y1  x  x_dot  phi  phi_dot  
Q_lqr = diag([ 1   1    1     1      1 ]);  % peso das variaveis 
R_lqr = 1;                                  % peso do esforço de controle

% Teste1
% Q_lqr = diag([0.1 0.1 100 0.1 0.5]); % peso alto para phi
% R_lqr = 1; % controle ainda penalizado normalmente
 
% Teste2
% Q_lqr = diag([1 1 100 1 1]);
% R_lqr = 10;

% Teste3
% Q_lqr = diag([1 10 100 10 1]);
% R_lqr = 1;

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
% figure; plot(traceP); ylabel('Traço de P - LQR'); xlabel('Amostras');

% Calculo do Ganho
K = ( Aa'*P*Ba*inv(Ba'*P*Ba+R_lqr) )';
disp('LQG K = '); disp(K);

% Analise dos poles e autovalores do sistema
disp('Autovalores do controlador LQR:');
disp(eig(Aa-Ba*K)); % precisam ser estaveis (está dentro do circulo únitario no plano z)

%% Estimador de estados - Filtro de Kalman
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
% figure; plot(tracePo); ylabel('Traço de P - Filtros de Kalman'); xlabel('Amostras');

% Calculo do Ganho
L = ( Aa*Po*Ca'*inv(Ca*Po*Ca'+R_kf) );
disp('KF L = '); disp(L);

% Analise dos poles e autovalores do sistema
disp('Autovalores do filtros de Kalman:');
disp(eig(Aa-L*Ca)); % precisam ser estaveis (está dentro do circulo únitario no plano z)

%% Simulação do Compensador dinâmico completo - LQG
t_f = 20;           % temos de simulação
N = round(t_f/Ts);  % número de amostras
t = 0:Ts:N*Ts-Ts;   % vetor de tempo discreto
    
    % Sinal de referencia
    r1(1:5) = 0; r1(6:N) = 0;

    % Disturbios na entrada e saida
    w1 = 0*wgn(1, N, 1e-3, 'linear'); % ruido
    w2(1:round(N/2)) = 0; w2(round(N/2)+1:N) = 1*0.0349; % deformação de ~2 graus sobre o pendulo

    v1 = 1*wgn(1, N, 1e-5, 'linear'); % ruido no sensor
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
    x(:,:, k) = Ad*x(:,:, k-1) + Bd*u(k-1) + [0; 0; 1; 0]*w2(k-1); 
    y(:,:, k) = Cd*x(:,:, k) + v1(k) + v2(k);

    % Modelo em espaço de estados aumentado com o filtro de Kalman
    xa(:,:, k) = Aa*xa(:,:, k-1) + Ba*du(k-1) + L*(y(:,:, k-1) - ya(:,:, k-1));
    ya(:,:, k) = Ca*xa(:,:, k);

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


%% Análise da estabilidade relativa - LQR
disp("Análise da estabilidade relativa - LQR");
% Método 1 - Análise em MA (recomendado apenas em sist. estáveis)
clear Tsen Ssen mt ms GmdB Pmdeg

figure; margin(ss(Aa, Ba, K, Da, Ts)); grid on;
title('Análise da estabilidade relativa em MA - LQR');

% Método 2 - Análise em MF (via funções de sensibilidade)
clear Tsen Ssen mt ms GmdB Pmdeg
    % Supondo o caso SISO, em que apenas r1 em r é um sinal não nulo, temos:
        % u = k1r1 - Kx
        % x_dot = Ax + Bu = Ax + B(k1r1 - Kx) = (A-BK)x + (Bk1)r1
        % y = Cx

Tsen = ss(Aa-Ba*K, Ba*K(1), Ca, Da, Ts);   % função sensibilidade complementar
Ssen = eye(ny,ny) - Tsen;                  % função sensibilidade

mt = max( sigma(Tsen) ); ms = max( sigma(Ssen) );

GmdB = min( 20*log10(ms/(ms-1)), 20*log10(1+(1/mt)) );
Pmdeg = (180/pi)*min( (2*asin(1/(2*ms)) ), (2*asin(1/(2*mt)) ) );

disp('GmdB = '); disp(GmdB); disp('Pmdeg = '); disp(Pmdeg);

figure; sigma(Tsen); hold; sigma(Ssen); grid;
legend('|Tsen|','|Ssen|'); title('Análise da estabilidade relativa em MF- LQR');


%% Análise da estabilidade relativa - Filtro de Kalman
disp("Análise da estabilidade relativa- KF");
clear Tsen Ssen mt ms GmdB Pmdeg
% Método 1 - Análise em MA (recomendado apenas em sist. estáveis)
    % x_dot_point = Ax_dot + Bu + Ly
    % y_dot = Cx_dot

figure; margin(ss(Aa, L, Ca, Da,Ts)); grid on;
title('Análise da estabilidade relativa em MA - KF');

% Método 2 - Análise em MF (via funções de sensibilidade)
    % x_dot_point = Ax_dot + Bu + L(y - y_dot) = (A-CL)x_dot + Bu + Ly
    % y_dot = Cx_dot
clear Tsen Ssen mt ms GmdB Pmdeg

Tsen = ss(Aa-L*Ca, L, Ca, Da, Ts);     % função sensibilidade complementar
Ssen = eye(ny,ny) - Tsen;              % função sensibilidade

mt = max( sigma(Tsen) ); ms = max( sigma(Ssen) );

% A margem é maior ou igual a isso, logo esse é o pior caso
GmdB = min( 20*log10(ms/(ms-1)), 20*log10(1+(1/mt)) );
Pmdeg = (180/pi)*min( (2*asin(1/(2*ms)) ), (2*asin(1/(2*mt)) ) );

disp('GmdB = '); disp(GmdB); disp('Pmdeg = '); disp(Pmdeg);

figure; sigma(Tsen); hold; sigma(Ssen); grid;
legend('|Tsen|','|Ssen|'); title('Análise da estabilidade relativa em MF- KF');


%% Análise da estabilidade relativa - LQG
disp("Análise da estabilidade relativa - LQG");
% Método 1 - Análise em MA (recomendado apenas em sist. estáveis)
clear Tsen Ssen mt ms GmdB Pmdeg A_comp B_comp C_comp

A_comp = [ Aa    zeros(n, n); 
          L*Ca   Aa-L*Ca];

B_comp = [Ba*K(1); Ba*K(1)];

C_comp = [Ca zeros(ny, n)];

Tsen = ss(A_comp, B_comp, C_comp, Da, Ts);   % função sensibilidade complementar
Ssen = eye(ny,ny) - Tsen;                    % função sensibilidade      

figure; sigma(Tsen); hold; sigma(Ssen); grid;
legend('|Tsen|','|Ssen|'); title('Análise da estabilidade relativa em MA - LQG');


% Método 2 - Análise em MF (via funções de sensibilidade)
clear Tsen Ssen mt ms GmdB Pmdeg A_comp B_comp C_comp

A_comp = [ Aa    -Ba*K; 
          L*Ca   Aa-Ba*K-L*Ca];

B_comp = [Ba*K(1); Ba*K(1)];

C_comp = [Ca zeros(ny, n)];

Tsen = ss(A_comp, B_comp, C_comp, Da, Ts);   % função sensibilidade complementar
Ssen = eye(ny,ny) - Tsen;                    % função sensibilidade      

mt = max( sigma(Tsen) ); ms = max( sigma(Ssen) );

GmdB = min( 20*log10(ms/(ms-1)), 20*log10(1+(1/mt)) );
Pmdeg = (180/pi)*min( (2*asin(1/(2*ms)) ), (2*asin(1/(2*mt)) ) );

disp('GmdB = '); disp(GmdB); disp('Pmdeg = '); disp(Pmdeg);

figure; sigma(Tsen); hold; sigma(Ssen); grid;
legend('|Tsen|','|Ssen|'); title('Análise da estabilidade relativa em MF- LQG');
     
