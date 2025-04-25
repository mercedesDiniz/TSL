%%%%
clear all; close all; clc;

%% Generating an input sequence to excite the underdamped process
tfinal = 4; % final time for the experiment, given in seconds
Ts = 0.05; % sampling time in seconds
N = round( tfinal/Ts ); % total number of samples