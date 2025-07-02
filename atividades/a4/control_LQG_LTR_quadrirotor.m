%% Atividade 4: Controle LQG/LTR para a velocidade lateral e longitudinal de um quadrirotor.
clear all; close all; clc;

%% Planta do sistema do tipo quadrirotor com saturação nos atuadores 
% Ref.: https://lacos.ufpa.br/plantas/ardrone_sim/ardrone_sim.html

path_lasse = 'C:\Users\Mercedes Diniz\Documents\PPGEE-UFPA\';
path_tsl = 'TSL-2025.1[Silveira]\atividades\a4\modelo\quadrotor_model.mat';

load([path_lasse, path_tsl]);

A = quadrotor_model.A 
B = quadrotor_model.B
C = quadrotor_model.C
D = quadrotor_model.D

Ts = quadrotor_model.Ts
