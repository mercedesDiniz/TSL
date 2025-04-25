clear all; close all; clc;

%% Experimental setup
tfinal = 20; % final time for the experiment, given in seconds
Ts = 0.05; % sampling time in seconds
N = round( tfinal/Ts ); % total number of samples


%% Loading the identified model
load('sys_model.mat');
Az = model.Az; 
    a1 = Az(2); a2 = Az(3);
Bz = model.Bz;
    b0 = Bz(2); b1 = Bz(3);


%% PID tuning
kp = 0.2; ki = 2; kd = 0.0;

    % Digital PID based on the Backward difference approximation
    s0 = kp +ki*Ts +kd/Ts;
    s1 = -kp -2*kd/Ts;
    s2 = kd/Ts;
    
    
    % Model-based PID tuning
    tau_mf = 0.5;
    s0 = ( 1-exp(-Ts/tau_mf) )/(b0+b1);
    s1 = a1*s0;
    s2 = a2*s0;
    
%% Control system analysis in the frequency domain
Gz = tf(Bz,Az,Ts);
Cz = tf([s0 s1 s2],[1  -1  0],Ts);
figure; margin(Cz*Gz); %title('Freq. Resp. of C(z)*G(z)');
Gmfz2 = feedback(Cz*Gz,1,-1);

%% Relative stability analysis in closed-loop
% The closed-loop system must be stable
disp('Closed-loop poles: '); pole(Gmfz2)

Tsen = Gmfz2; % co-sesitivity function
Ssen = 1 -Tsen;

mt = max( sigma(Tsen) );
ms = max( sigma(Ssen) );

GmdB = min( 20*log10(ms/(ms-1)), 20*log10(1+(1/mt)) );
Pmdeg = (180/pi)*min( (2*asin(1/(2*ms)) ), (2*asin(1/(2*mt)) ) );

disp('Gain and Phase margins obtained by closed-loop analysis:');
disp('GmdB = '); disp(GmdB);
disp('Pmdeg = '); disp(Pmdeg);

figure; sigma(Tsen); hold; sigma(Ssen); grid;
legend('|Tsen|','|Ssen|');

%% Closed-loop simulation

% Generating a reference sequence
r(1:4)=0; r(10:1*N/4)=1; r(1*N/4 +1: 2*N/4)=3; 
r(2*N/4 +1:3*N/4)=4; r(3*N/4 +1: 4*N/4)=2; % given in Volts

% Load or noise disturbances
  v(1:round(N/2))=0; v(1+round(N/2):N)=0;

daqduino_start('/dev/ttyUSB0');
for k = 1:2
    ym(k) = 0; um(k) = 0; em(k) = 0;
    y(k) = 0; u(k) = 0; e(k) = 0;
end
for k = 3:N
    tic
    % Simulation part
    ym(k) = -a1*ym(k-1) -a2*ym(k-2) +b0*um(k-1) +b1*um(k-2) ...
        +v(k)+a1*v(k-1)+a2*v(k-2);
    
        % PID control
        em(k) = r(k)-ym(k);
        um(k) = um(k-1) +s0*em(k) +s1*em(k-1) +s2*em(k-2);
        
        % Control saturation
        if um(k) <= 0
            um(k) = 0;
        elseif um(k) >= 5
            um(k) = 5;
        end
    
    %% Real control system
    y(k) = daqduino_read();
        % Real PID control
        e(k) = r(k)-y(k);
        u(k) = u(k-1) +s0*e(k) +s1*e(k-1) +s2*e(k-2); 
    daqduino_write( u(k), Ts );
    toc
end
daqduino_end;

%% Plots
t = 0:Ts:N*Ts-Ts;
figure;
subplot(211)
    plot(t,r,'k',t,ym,'b',t,y,'r');
    legend('Ref.','Sim','Real');
ylabel('Amplitude (V)'); xlabel('Time (s)');
subplot(212)
    plot(t,um,'b',t,u,'r');
    legend('Sim','Real');
ylabel('Control (V)'); xlabel('Time (s)');