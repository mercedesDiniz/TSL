% Exemplo: linearizacao local
% Obter uma aproximacao de y(u) = sqrt(u)
clear all; close all; clc;

% Considerando a funcao nao linear
u = 0:0.01:60;
N = length(u); % Numero de iteracoes pra montar o vetor

% y(k), k=1,2,...,N, equivalente a y(u).
for k=1:N
y(k) = sqrt(u(k));
end
plot(u,y); ylabel('y(u)'); xlabel('u'); hold;

% Aproximacao linear: L(u) = f(ueq) +slope*(u - ueq)
% Selecao do ponto em u onde se quer linearizar a funcao
ueq = 40; % vamos obter uma aproximacao linear L(u)
% em torno de ueq.
% Apoio ao calculo diferencial usando Matlab:
% via toolbox symbolic e funcao diff:
syms usym; % gera uma variavel simbolica usym para determinar a
% derivada da funcao y(usym) = sqrt(usym)
disp('Inclinacao da reta');
slope = diff( sqrt(usym) )
pretty(slope); % apresenta de uma forma legivel

% Reultado do slope = 1/( 2*sqrt(ueq) )
for k=1:N
% L(u) = f(ueq) +slope*(ueq -x)
L(k) = sqrt(ueq) +( 1/(2*sqrt(ueq)) )*( u(k) - ueq );
end
plot(u,L);
legend('Nao linear','Aproximacao','Location','NorthWest');
title(['Aproximacao em torno de ueq = ' num2str(ueq)]);