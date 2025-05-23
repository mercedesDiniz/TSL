%% Problema exemplo de controle via realimentação de estado
clear all; close all; clc;

wn = 2;
zeta = 0.3;

A = [0 1;  -wn^2  -2*zeta*wn]
B = [0; wn^2]
C = [1 0]
D = [0]
sys = ss(A,B,C,D)
