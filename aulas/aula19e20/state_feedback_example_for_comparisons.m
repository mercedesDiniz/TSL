%% Servo control of a 2nd order system based on
% state feedback control.
% The control system must guarantee a null offset error
% and step-like load disturbances rejection.
clear all; close all; clc;

%% System model in the continuous freq. domain
ks = 2; % V/V
wn = 10; % rad/s
zeta = 0.1; % damping factor

Gs = tf([ks*wn^2],[1    2*zeta*wn  wn^2])

%% Converting the system to a state space realization
Ac = [0  1;   -wn^2  -2*zeta*wn];
Bc = [0;  ks*wn^2];
C = [1 0];
D = [0];
sysc = ss(Ac,Bc,C,D);

%% Converting the stste space system to the discrete
% time domain k:
Ts = 0.05; % sampling time in seconds
sysd = c2d(sysc,Ts);
[A,B,C,D] = ssdata(sysd) % discrete model matrices

%% Augmenting the discrete time state space model
% to include the discrete integrator.
Aa = [  1     C*A;
      [0;0]    A  ];
Ba = [C*B; B];
Ca = [1 0 0];
Da = 0;
sysa = ss(Aa,Ba,Ca,Da,Ts);

%% Control system design using Ackermann's formula
Co = [Ba Aa*Ba (Aa^2)*Ba]; % Co = ctrb(Aa,Ba)

    % If the det(Co) is different from zero, then the
    % realization is Controllable.
    % For MIMO systems, the controlability check must be
    % made using the rank(Co), which must be equal to
    % the order of the system.
    if rank(Co) == 3
        disp('The system is controllable');
    else
        disp('The system is not controllable');
    end
    
    % Specifying the closed-loop dynamics
    w_cl = wn/2; % rad/s
    zeta_cl = 1; % closed-loop damping
        w_slack = 100; % slack frequency
    Gcl = tf(1,conv([1   2*zeta_cl*w_cl  w_cl^2], ...
        [1 w_slack]));
    Gcld = c2d(Gcl,Ts);
       Pcl = Gcld.den{1};
    
    % Ackermann's formula
    K = [0 0 1]*inv(Co)*(Aa^3 +Pcl(2)*Aa^2 ...
                        +Pcl(3)*Aa +Pcl(4)*Aa^0);
    K = acker(Aa,Ba,roots(Pcl) );
    
    
%% Observer system design using Ackermann's formula
Ob = [Ca; Ca*Aa; Ca*(Aa^2)]
Ob = obsv(Aa,Ca)

    % If the det(Ob) is different from zero, then the
    % realization is Observable.
    % For MIMO systems, the observability check must be
    % made using the rank(Ob), which must be equal to
    % the order of the system.
    if rank(Ob) == 3
        disp('The system is Observable');
    else
        disp('The system is not Observable');
    end
    
    % Specifying the closed-loop dynamics of the OBSERVER
    w_clo = wn/2; % rad/s
    zeta_clo = 1; % closed-loop damping
        w_slacko = 100; % slack frequency
    Gclo = tf(1,conv([1   2*zeta_clo*w_clo  w_clo^2], ...
        [1 w_slacko]));
    Gcldo = c2d(Gclo,Ts);
       Pclo = Gcldo.den{1};
    
    % Ackermann's formula for the observer
    L = [0 0 1]*inv(Ob')*(Aa'^3 +Pclo(2)*Aa'^2 ...
                        +Pclo(3)*Aa' +Pclo(4)*Aa'^0); L=L'
    L = (acker(Aa',Ca',roots(Pclo) ))'
    
%% Control system simulation
tfinal = 10; % total simulation time in seconds
N = round( tfinal/Ts ); % total number of samples
t = 0:Ts:N*Ts-Ts; % time vector for the plots

    r1(1:5)=0; r1(6:N)=0;

    % Disturbances at the output
      w = 0*wgn(1,N, 1e-3, 'linear'); % Noise
      w2(1:round(N/2))=0; w2(round(N/2)+1:N)=+0.02;
    
    % Disturbances at the output
      v = 1*wgn(1,N, 1e-3, 'linear'); % Noise
      v2(1:round(N/2))=0; v2(round(N/2)+1:N)=+0.2*0;

for k = 1:1
    % Nominal model initial conditions
    x(:,:,k) = [0.1;0];
    y(k) = C*x(:,:,k);
    u(k) = 0;
    
    % Augmented model initial conditions
    xa(:,:,k) = [0;0;0];
    ya(k) = Ca*xa(:,:,k);
    du(k) = 0;
    
        % Estimated state variables based on the augmented ones
        x1e(k)=0;
        x2e(k)=0;
end
for k = 2:N
    % State space model simulation
    x(:,:,k) = A*x(:,:,k-1) +B*u(k-1) +[1;0]*w2(k-1);
    y(k) = C*x(:,:,k) +v(k) +v2(k);
    
    
    % Augmented state space model
    xa(:,:,k) = Aa*xa(:,:,k-1) +Ba*du(k-1) +L*( y(k-1)-ya(k-1) ) ;
    ya(k) = Ca*xa(:,:,k);
        % Estimated state variables based on the augmented ones
            x1e(k)=x1e(k-1) +xa(2,k);
            x2e(k)=x2e(k-1) +xa(3,k);
    
    % State-feedback control law
    du(k) = K*( [r1(k);0;0] -xa(:,:,k) );
    
    % Passing du(k) through the discrete integrator
    u(k) = u(k-1) +du(k); %u(k)=r1(k);
end
    

figure; % I/O plots
subplot(211)
    plot(t,r1,'k'); hold;
    plot(t,y,'b'); plot(t,ya,':r');
    legend('r1(t)','y(t)','ya(t): estimated');

subplot(212)
    plot(t,u,'b');
    
figure; % estimated state variables of the augmented estimator
subplot(311)
    plot(t,xa(1,:),'r');
    ylabel('xa_1(t)');
    title('State variables of the augmented estimator');

subplot(312)
    plot(t,xa(2,:),'r');
    ylabel('xa_2(t)');

subplot(313)
    plot(t,xa(3,:),'r');
    ylabel('xa_3(t)');

    
figure; % estimated state variables based on the nominal system
subplot(211)
    plot(t,x(1,:),'b'); hold;
    plot(t,x1e,'r');
    ylabel('xa_1(t)');
    title('Estimated state variables based on the nominal system');

subplot(212)
    plot(t,x(2,:),'b'); hold;
    plot(t,x2e,'r');
    ylabel('xa_2(t)');




