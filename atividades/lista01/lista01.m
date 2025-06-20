%% Lista 01
clear all; close all; clc;

%% 1) Calcule os determinantes das matrizes:
disp("QUESTÃO 1");
clear all; 
M1 = [1 3 0; -1 2 1; 1 -1 1] 
disp("det(M1) = "); det(M1)

M2 = [-1 3 2; -1 2 1; 2 -6 -4]
disp("det(M2) = "); det(M2)

M3 = [1 0 3 5; 3 1 1 0; 0 1 1 1; 1 2 3 4]
disp("det(M3) = "); det(M3)

M5 = [-1 3 0 0 0 0; 2 1 0 0 0 0; 0 0 1 2 0 0; 0 0 -1 -1 0 0; 0 0 0 0 5 6; 0 0 0 0 1 1]
disp("det(M5) = "); det(M5)

%% 2) Quando possível, calcule as inversas das seguintes matrizes:
disp("QUESTÃO 2");
clear all;
A = [1 2; 3 1]
disp("inv(A) = "); inv(A)

B  = [ 1 3 0; -1 2 1; 1 -1 1]
disp("inv(B) = "); inv(B)

C = [-1 3 2; -1 2 1; 2 -6 -4]
disp("inv(C) = Não Existe. "); inv(C) 

%% 8) Calcule os autovalores das seguintes matrizes:
disp("QUESTÃO 8");
clear all; 

A = [1 2; 0 4]
disp("eig(A) = "); eig(A)

B = [1 2; 2 4]
disp("eig(B) = "); eig(B)

C = [2 2; 8 8]
disp("eig(C) = "); eig(C)

D = [1 0; 1 4]
disp("eig(D) = "); eig(D)

%%  9) Trace os Diagramas de Bode dos seguintes sistemas:
disp("QUESTÃO 9");
clear all; 

% Letra (a)
G_a = tf(1, 1);

% Letra (b)
G_b = tf(1, [1 0]);

% Letra (c)
G_c = tf(1, [1 1]);

% Letra (d)
G_d = tf(4, conv([1 3], [1 8])); % (s+3)(s+8)

% Letra (e)
num_e = 5 * conv([0.01 1], [0.01 1]);               % (1 + 0.01s)^2
den_e = conv([1 0], conv([0.0001 1], [0.0001 1]));  % s*(1 + 0.0001s)^2
G_e = tf(num_e, den_e);

% Letra (f)
G_f = tf(4, conv([1 5], [1 5])); % (s+5)^2

% Plots
figure;
subplot(3,2,1);
bode(G_a); grid on; title('G_a(s) = 1');

subplot(3,2,2);
bode(G_b); grid on; title('G_b(s) = 1/s');

subplot(3,2,3);
bode(G_c); grid on; title('G_c(s) = 1/(s+1)');

subplot(3,2,4);
bode(G_d); grid on; title('G_d(s) = 4/[(s+3)(s+8)]');

subplot(3,2,5);
bode(G_e); grid on; title('G_e(s) = 5(1+0.01s)^2 / [s(1+0.0001s)^2]');

subplot(3,2,6);
bode(G_f); grid on; title('G_f(s) = 4/(s+5)^2');

%% 13)
% Elabore um exemplo com um processo sub-amortecido de segunda ordem 
% ligado à uma malha de controle PI. Usando este exemplo, avalie as margens
% de Ganho e de Fase pelo método de Bode (em malha aberta) e compare com a
% avaliação a partir das margens obtidas via funções de sensibilidade 
% (análise em malha fechada).

disp("QUESTÃO 13");
clear all; 

% 1. Planta subamortecida de 2ª ordem
wn = 4;        % frequência natural
zeta = 0.3;    % fator de amortecimento
G = tf(wn^2, [1 2*zeta*wn wn^2]);

% 2. Controlador PI
Kp = 1.2; Ki = 1.5;
C = tf([Kp Ki], [1 0]);  % C(s) = (Kp*s + Ki)/s

% 3. Análise em malha aberta (C(s)*G(s))
L = C * G;  % Malha aberta

figure;
margin(L); grid on;
title('Diagrama de Bode - Análise em Malha Aberta');

% 4. Malha fechada com realimentação negativa unária
Tsen = feedback(L, 1);      % Função de co-sensibilidade T(s)
Ssen = 1 - Tsen;            % Função de sensibilidade S(s)

mt = max( sigma(Tsen) );
ms = max( sigma(Ssen) );

GmdB = min( 20*log10(ms/(ms-1)), 20*log10(1+(1/mt)) );
Pmdeg = (180/pi)*min( (2*asin(1/(2*ms)) ), (2*asin(1/(2*mt)) ) );

disp('--- Análise com funções de sensibilidade ---');
disp(['Ganho máximo |Tsen| (mt): ', num2str(mt)]);
disp(['Ganho máximo |Ssen| (ms): ', num2str(ms)]);
disp(['Margem de ganho estimada (dB): ', num2str(GmdB)]);
disp(['Margem de fase estimada (graus): ', num2str(Pmdeg)]);

figure; sigma(Tsen); hold; sigma(Ssen); grid;
legend('|Tsen|','|Ssen|'); title('Funções de sensibilidade - Análise em Malha Fechada');

%% 14)
% Utilizando o mesmo exemplo elaborado na Questão 13, obtenha uma
% realização controlável em espaço de estados para este sistema e projete 
% um regulador por realimentaçãode estados tal que os polos de malha 
% fechada sejam puramente reais. 
% Após, mostre como é possível verificar o Diagrama de Bode da malha direta
% utilizando as matrizes do sistema emespaço de estados e o ganho do regulador.

disp("QUESTÃO 14");

% 1. Realização controlável
A = [0 1; -wn^2 -2*zeta*wn];
B = [0; 1];
C = [wn^2 0];
D = 0;

sys = ss(A,B,C,D)

% 2. Projeto de realimentação de estados
% Desejamos polos puramente reais
p_des = [-4 -6];  % exemplo

K = place(A, B, p_des)

% Sistema em malha fechada com realimentação de estados (A - B*K)
Acl = A - B*K;

% Verificar estabilidade
disp('Polos de malha fechada com realimentação de estados:');
disp(eig(Acl));

% 3. Malha direta L(s) = K*(sI - A)^(-1)*B (sem realimentação unária)
[num, den] = ss2tf(A, B, K, 0);

L = tf(num, den)

% 4. Diagrama de Bode da malha direta
figure; margin(L); grid on;
title('Diagrama de Bode da malha direta K*(sI - A)^(-1)*B');

%% 17)
% Para o sistema mostrado a seguir, utilizando a fórmula de Ackermann, 
% projete um observador de estados tal que os pólos do sistema observador 
% sejam o dobro dos pólos de malha aberta da planta. 
% Apresente, como resultado, o ganho L do observador.

disp("QUESTÃO 17");
clear all;

% Matrizes
A = [1 2; -3 -4];
B = [1; 2];
C = [1 0];

% 1. Autovalores da planta
p_plant = eig(A)

% 2. Dobro dos polos da planta 
p_observer = 2 * p_plant;

% 3. Ganho do observador (Ackermann)
disp('Ganho do observador L:')
L = acker(A', C', p_observer)'

%% 19)
% 19) Dado o sistema mostrado a seguir, descrito no espaço de estados, 
% verifique se este sistema é (i) estável, se é (ii) controlável e (iii)
% observável. Por fim, (iv) obtenha a função detransferência que descreve 
% a sua relação de entrada e saída.
disp("QUESTÃO 19");
clear all;

% Matrizes
A = [20 -200; 0 0];
B = [0; 1];
C = [0 200];
D = 0;

% Verificar se o sistema é (i) estável:
autovalores = eig(A);
if all(real(autovalores) < 0)
    disp('Sistema Estável');
else
    disp('Sistema NÃO é Estável');
end

% Verificar se o sistema é (ii) controlável:
Co = ctrb(A, B);
if rank(Co) == size(A,1)
    disp("Sistema Controlável");
else
    disp("Sistema NÃO é Controlável");
end

% Verificar se o sistema é (iii) observável:
Ob = obsv(A, C);

if rank(Ob) == size(A,1)
    disp("Sistema Observável");
else
    disp("Sistema NÃO é Observável");
end

% Calcular a função de transferência
s = tf('s');  % Define a variável 's' para a função de transferência
G = C * inv(s * eye(size(A)) - A) * B + D

%% 21
% Dado o sistema mostrado a seguir, descrito no espaço de estados, 
% verifique se este sistemaé (i) estável, se é (ii) controlável e (iii) 
% observável. Por fim, (iv) obtenha a função detransferência que descreve 
% a sua relação de entrada e saída.

disp("QUESTÃO 21");
clear all;

% Matrizes
A = [3 -2 0; 1 0 0; 0 1 0];
B = [1; 0; 0];
C = [0 0 1];
D = 0;

% Verificar se o sistema é (i) estável:
autovalores = eig(A);
if all(real(autovalores) < 0)
    disp('Sistema Estável');
else
    disp('Sistema NÃO é Estável');
end

% Verificar se o sistema é (ii) controlável:
Co = ctrb(A, B);
if rank(Co) == size(A,1)
    disp("Sistema Controlável");
else
    disp("Sistema NÃO é Controlável");
end

% Verificar se o sistema é (iii) observável:
Ob = obsv(A, C);

if rank(Ob) == size(A,1)
    disp("Sistema Observável");
else
    disp("Sistema NÃO é Observável");
end

% Calcular a função de transferência
s = tf('s');  % Define a variável 's' para a função de transferência
G = C * inv(s * eye(size(A)) - A) * B + D