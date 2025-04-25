%% System Identification
clear all; close all; clc;

%% Experimental setup
tfinal = 10; % final time for the experiment, given in seconds
Ts = 0.05; % sampling time in seconds
N = round( tfinal/Ts ); % total number of samples

%% Generating an input sequence to excite the underdamped process 
u(1:4)=0; u(10:1*N/4)=1; u(1*N/4 +1: 2*N/4)=3; u(2*N/4 +1:3*N/4)=4;
u(3*N/4 +1: 4*N/4)=2; % given in Volts

daqduino_start('/dev/ttyUSB0');
for k = 1:N
    tic
    y(k) = daqduino_read();
    
    daqduino_write( u(k), Ts );
    toc
end
daqduino_end;

%% Least Squares Parametric Estimation
tini = 2.5;
Nini = round(tini/Ts);
tfin = 3.5;
Nfin = round(tfin/Ts);

ui=u(Nini:Nfin);
yi=y(Nini:Nfin);
Ni = length(yi);

PHI=[];
for k=3:Ni
    PHI = [PHI ; [-yi(k-1) -yi(k-2) ui(k-1) ui(k-2)] ];
end
theta = inv(PHI'*PHI)*PHI'*yi(3:Ni)'
a1 = theta(1); a2 = theta(2); b0 = theta(3); b1 = theta(4);

ym(1:2)=y(1:2);
for k = 3:N
    ym(k) = -a1*ym(k-1) -a2*ym(k-2) +b0*u(k-1) +b1*u(k-2);
end

%% Plots
t = 0:Ts:N*Ts-Ts;
plot(t,u,'r',t,y,'b',t,ym,'k'); legend('u(t)','y(t)','ym(t)');
ylabel('Amplitude (V)'); xlabel('Time (s)');

%% Saving ident data and model data
% save 'datalog_tuy.txt' datalog -ascii

model.Az = [1 a1 a2];
model.Bz = [0 b0 b1];
model.Ts = Ts;
save 'sys_model.mat' model;


