%% Atividade 02: Controle do pêndulo invertido sobre carro usando LQG.
clear all; close all; clc;

%% Planta do pêndulo invertido 
% Ref.: https://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=SystemModeling 
% A variável sensorialmente medida é a posição angular do pêndulo (phi)

% Parametros do sistema
M = 0.5; m = 0.2; b = 0.1; I = 0.006; g = 9.8; l = 0.3;   

% Modelo em Espaço de Estados Contínuo
p = I*(M+m)+M*m*l^2; 

A = [0      1              0           0;
     0 -(I+m*l^2)*b/p  (m^2*g*l^2)/p   0;
     0      0              0           1;
     0 -(m*l*b)/p       m*g*l*(M+m)/p  0];

B = [   0; 
    (I+m*l^2)/p;
        0; 
      m*l/p];

C = [0 0 1 0]; % só mede phi
     
D = 0;

sys_ss_c = ss(A,B,C,D, ...
    'statename', {'x' 'x_dot' 'phi' 'phi_dot'}, ...
    'inputname',{'u'}, ...
    'outputname', {'phi'})

% Modelo de Espaço de Estado Discreto
Ts = 0.01; % s
sys_ss_d = c2d(sys_ss_c, Ts)
[Ad, Bd, Cd, Dd] = ssdata(sys_ss_d);

% Modelo Aumentando de Espaço de Estado Discreto
% Adicionar o integrador
Aa = [1,            Cd*Ad;
      zeros(4,1),   Ad];  

Ba = [Cd*Bd;
       Bd];

Ca = [Cd, 0];

Da = Dd;

sys_ss_d_a = ss(Aa, Ba, Ca, Da, Ts,  ...
    'statename', {'x' 'x_dot' 'phi' 'phi_dot' 'int'}, ...
    'inputname',{'u'}, ...
    'outputname', {'phi'})


% Verificação de controlabilidade e observabilidade
Co = ctrb(Ad, Bd); Ob = obsv(Ad, Cd);

if rank(Co) == size(Ad,1)
    disp("Sistema Controlável");
else
    warning("Sistema NÃO é Controlável");
end

if rank(Ob) == size(Ad,1)
    disp("Sistema Observável");
else
    warning("Sistema NÃO é Observável");
end

%% Controlador LQR
             % x  x_dot  phi  phi_dot  int     % inclui o integrador
Q_lqr = diag([ 1    1    10    1     1]);      % peso das variaveis de estados
R_lqr = 1;                                     % peso do esforço de controle

K = dlqr(Aa, Ba, Q_lqr, R_lqr)

% %% Analise da resposta em frequencia do LQR
% Tsen = ss(Aa-Ba*K,Ba*K(1),Ca,Da,Ts); % TODO: K(3) ou K(1)?
% Ssen = eye(1,1)-Tsen; 
% 
% mt = max( sigma(Tsen) );
% ms = max( sigma(Ssen) );
% 
% GmdB = min( 20*log10(ms/(ms-1)), 20*log10(1+(1/mt)) );
% Pmdeg = (180/pi)*min( (2*asin(1/(2*ms)) ), (2*asin(1/(2*mt)) ) );
% 
% disp('Análise LQG');
% disp('Margens de ganho e fase obtidas por análise de malha fechada:');
% disp('GmdB = '); disp(GmdB);
% disp('Pmdeg = '); disp(Pmdeg);
% 
% figure; sigma(Tsen); hold; sigma(Ssen); grid;
% legend('|Tsen|','|Ssen|'); title('LQR Sensitivities');

%% Filtro de Kalman (Observador)
            % x  x_dot  phi  phi_dot  % estimando apenas os estados físicos
Q_kf = diag([ 1    1    10     1 ]);  % pesos das variaveis de estados
R_kf = 1;                             % ruído da medição  
L = (dlqr(Ad',Cd',Q_kf,R_kf))'        % ganho do observador 

% %% Analise da resposta em frequencia do Filtro de Kalman
% % clear Tsen Ssen mt ms GmdB Pmdeg
% Tsen = ss(Aa-L*Ca,L,Ca,Da,Ts);
% Ssen = eye(1,1)-Tsen;
% 
% mt = max( sigma(Tsen) );
% ms = max( sigma(Ssen) );
% 
% GmdB = min( 20*log10(ms/(ms-1)), 20*log10(1+(1/mt)) );
% Pmdeg = (180/pi)*min( (2*asin(1/(2*ms)) ), (2*asin(1/(2*mt)) ) );
% 
% disp('Análise do Filtro de Kalman');
% disp('Margens de ganho e fase obtidas por análise de malha fechada:');
% disp('GmdB = '); disp(GmdB);
% disp('Pmdeg = '); disp(Pmdeg);
% 
% figure; sigma(Tsen); hold; sigma(Ssen); grid;
% legend('|Tsen|','|Ssen|'); title('Kalman Filter Sensitivities');
% 
% %% Compensador dinâmico completo LQG
% 
% % Configurações da simulação
% tfinal = 10;                % tempo de simulação (s)
% N = round(tfinal/Ts);       % numero de amostras
% t = 0:Ts:N*Ts-Ts;           
% 
% % Condições iniciais
% x(:,1) = [0.1; 0; 0.1; 0];   % posição, velocidade, angulo, velocidade angular
% xa(:,1) = [0; 0; 0; 0; 0];   % integrador + estados estimados
% 
% y(1) = C*x(:,1);             % saida medida inicial
% ya(1) = Ca*xa(:,1);          % saida estimada inicial
% 
% x1e(1) = xa(2,1); 
% x2e(1) = xa(3,1);
% x3e(1) = xa(4,1);
% x4e(1) = xa(5,1);
% 
% % Disturbios na entrada e saida
% w = 1*wgn(1,N, 1e-3, 'linear');
% % w2(1:round(N/2))=0; w2(round(N/2)+1:N)=+0.02; 
% v = 1*wgn(1,N, 1e-3, 'linear'); 
% v2(1:round(N/2))=0; v2(round(N/2)+1:N)=+0.2*0;
% 
% ref(1:round(N/6))= 0; ref(round(N/6):N)=1;   % referencia 
% u = zeros(1, N);    % controle aplicado
% du = zeros(1, N);  
% 
% % Simulação
% for k = 2:N
%     % Sistema real
%     x(:,k) = Ad*x(:,k-1) + Bd*u(k-1) + [1;0;0;0]*w(k); 
%     y(k) = C*x(:,k) + v(k) +v2(k);                        
%     
%     % Estimador de estados (Filtro de Kalman)
%     xa(:,k) = Aa*xa(:,k-1) + Ba*du(k-1) + L*( y(k-1) - ya(k-1) );
%     ya(k) = Ca*xa(:,k);
% 
%     % Controle LQR com ação integrativa
%     du(k) = K * ( [ref(N); zeros(4,1)] - xa(:,k) ); % controle com integrador
%     
%     % Atualiza controle total aplicando integrador na entrada
%     u(k) = u(k-1) + du(k);
% 
%     % Estados estimados individualmente
% %     x1e(k) = xa(2,k);
% %     x2e(k) = xa(3,k);
% %     x3e(k) = xa(4,k);
% %     x4e(k) = xa(5,k);
% end
% 
% % Calculo da Função Custo (TODO: validacao pendente)
% % J = sum(x(3,:).^2 + u.^2)
% 
% %% Plots
% % Comparação entre referência, saída real e estimada
% figure;
% subplot(2,1,1)
%     plot(t, ref, 'k--'); hold on;
%     plot(t, y, 'b'); 
%     plot(t, ya, 'r:');
%     legend('Ref', '\phi (real)', '\phi (estimada)');
%     ylabel('\phi (rad)');
%     title('Seguimento de referência e estimativa');
% 
% subplot(2,1,2)
%     plot(t, u, 'm');
%     ylabel('u(t) [N]');
%     xlabel('Tempo [s]');
%     title('Sinal de controle');
