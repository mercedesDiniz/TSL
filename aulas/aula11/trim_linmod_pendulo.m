%% Realiza a trimagem do modelo não linear do pêndulo invertido sobre o
%% carro para encontrar uma solução de equilíbrio para a etapa de
%% linearização assistida por computador.
clear all; close all; clc;


%% Busca solução via método iterativo de trimagem (trim)
x0 = [0 0 pi 0]'; %[xp d/dt_xp theta d/dt_theta]
F0 = 0;

[Xeq,Feq,Yeq,DXeq]=trim('inv_pend_R2014a',x0,F0)


%% Lineariza o modelo usando Xeq e Feq determinados
[A,B,C,D]=linmod('inv_pend_R2014a',Xeq,Feq)


%% Confirmando os valores do modelo linear obtido com aqueles
%% calculados manualmente:
M = 0.5;
m = 0.2;
b = 0.1;
l = 0.3;
I = 0.006;
d_en = (M+m)*(I+m*l*l)-(m*l)*(m*l);
g = 9.8;

Am = [ 0  1                  0                   0;
       0 -(b*(I+m*l^2))/d_en (((m*l)^2)*g)/d_en  0;
       0  0                  0                   1;
       0 -(m*l*b)/d_en       (m*l*g*(M+m))/d_en  0];
Bm = [0; (I+m*l^2)/d_en ; 0 ; m*l/d_en ];

A,Am
B,Bm


%% Análise do modelo linearizado
sys = ss(A,B,[1 0 0 0;0 0 1 0],0);
disp('Ganho Estático = '); dcgain(sys)
bode(sys);
disp('Autovalores = '); eig(A)

