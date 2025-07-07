%% LQG Control Applied to an Inverted Pendulum on a Cart
%% SIDO case: controlling the cart position while regulating the
% pendulum.
clear all; close all; clc;

%% Inverted Pendulumn linear model in continuous time domain
    % Parameters
    M = 0.5; % Cart mass in kg
    m = 0.2; % Pendulum mass in kg
    b = 0.1; % 
    l = 0.3;
    I = 0.006;
    d_en = (M+m)*(I+m*l*l)-(m*l)*(m*l);
    g = 9.8;

    % State equation matrices
    Am = [ 0  1                  0                   0;
           0 -(b*(I+m*l^2))/d_en (((m*l)^2)*g)/d_en  0;
           0  0                  0                   1;
           0 -(m*l*b)/d_en       (m*l*g*(M+m))/d_en  0];
    Bm = [0; (I+m*l^2)/d_en ; 0 ; m*l/d_en ];
    
    % Cart control input: u = F: Force in Newtons
    
    % State vector: x^T = [ xp dot_xp theta dot_theta ]
        % x1 = xp: cart position in meters
        % x2 = dot_xp: cart velocity in m/s
        % x3 = theta: Pendulum angle (rad)
        % x4 = dot_theta: Pendulum angular velocity (rad/s)

    % Selecting x1 and x3 as sensor-measured outputs (C matrix)
    C = [1 0 0 0;0 0 1 0]; ny = 2; % number of outputs
         % y1 = x1: cart position (m)
         % y2 = x3: Pendulum angle (rad)
         
    D=zeros(ny,1);
         
    % State space system in continuous time domain
    sysc = ss(Am,Bm,C,D);
    
        % System poles/eigenvalues
        disp('----------------------------------------------------------');
        disp('Unstable system with the following continuous eigenvalues:');
        disp(eig(Am));
        disp('----------------------------------------------------------');
        
        % Checking max frequencies to select the sampling frequency
        disp('----------------------------------------------------------');
        disp('Assessing max frequency of the system:');
        damp(sysc)
        fmax = max(damp(sysc));
        disp(['Maximum frequency fmax = ' num2str(fmax) ' rad/s']);
        disp('----------------------------------------------------------');
        
    % Discrete-time state space system
        vfactor = 200; % sampling velocity factor
        ws = vfactor*fmax;
            % Remark: rules for sampling open-loop stable systems may
            %         not apply to unstable systems cases.
        fs = ws/(2*pi);
        Ts = 1/fs % seconds
        Ts = fix(vfactor*Ts)/vfactor; % Rounding
        disp('----------------------------------------------------------');
        disp(['Sampling time Ts = ' num2str(Ts) ' seconds']);
        disp('----------------------------------------------------------');
        
        sysd = c2d(sysc,Ts);
        [A,B,C,D] = ssdata(sysd);
        
    % Comparing the continuous model to the discrete one in the freq. dom.
    figure; bode(sysc,sysd); title('SIDO Bode Plots');
    legend('Continuous','Discrete','Location','SouthWest');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% LQG Design
    % Augmenting the state space system to include integrators
        Aa = [  eye(ny,ny)     C*A;
               zeros(4,ny)    A  ];
        Ba = [C*B; B];
        Ca = [eye(ny,ny) zeros(ny,4)];
        Da = D;
        sysa = ss(Aa,Ba,Ca,Da,Ts);
        
        % Assessing the augmented system in the frequency domain
        figure; bode(sysa); grid; title('Augmented SIDO model');
        
	% Controllability and Observability check
    n = length(Aa); % system order
    Co = ctrb(Aa,Ba); % Controllability matrix
    Ob = obsv(Aa,Ca); % Observability matrix
    
    disp('----------------------------------------------------------');
    if rank(Co) == n
        disp('rank(Co) = n and the augmented realization is controllable');
    else
        disp('rank(Co) =/= n and the realization is not controllable');
    end
        if rank(Ob) == n
        disp('rank(Ob) = n and the augmented realization is observable');
    else
        disp('rank(Ob) =/= n and the realization is not observable');
    end
    disp('----------------------------------------------------------');
    
    disp('----------------------------------------------------------');
    disp(['REMARK: when the system is not controllable or not       ';
          'observable, please refer to the section                  ';
          '|Selection of the PI Weighting Matrices| in the book of  ';
          'Stevens and Lewis (2016), on page 405.                   ';
          'When this happens, we assess DETECTABILITY from the pair ';
          'of matrices sqrt(Q) and A.                               ']);
    disp('----------------------------------------------------------');
    
    
    % LQR Performance Index Weightings on Qlq and Rlq
                 % xa1 xa2  dx1 dx2 dx3 dx4
        Qlq = diag([1   1    1e3   1e3   1   1]);
        Rlq = 1e3;
    disp('----------------------------------------------------------');
    disp('Checking the DETECTABILITY of (sqrt(Qlq),Aa):');
    disp('Ob = obsv(Aa,sqrt(Qlq))');
    Ob_lq = obsv(Aa,sqrt(Qlq));
    
    % Kalman Filter Performance Index Weightings on Qkf and Rkf
                  % xa1  xa2    dx1 dx2 dx3 dx4
        Qkf = diag([1e3  1e3    1    1   1   1]);
        Rkf = diag([1e3 1e3]);
    disp('----------------------------------------------------------');
    disp('Checking the DETECTABILITY of (sqrt(Qkf)^T,Aa^T):');
    disp('Ob = obsv(Aa^T,sqrt(Qlq)^T)');
    Ob_kf = obsv(Aa',sqrt(Qlq)');
    
    disp('----------------------------------------------------------');
    if rank(Ob_lq) == n
        disp('rank(Ob_lq) = n, system is detectable from the LQR problem');
    else
        disp('rank(Ob_lq) =/= n, system IS NOT detectable from the LQR problem');
    end
        if rank(Ob_kf) == n
        disp('rank(Ob_kf) = n, system is detectable from the KF problem');
    else
        disp('rank(Ob_kf) =/= n, system IS NOT detectable from the KF problem');
    end
    disp('----------------------------------------------------------');
    
    
    %% LQR Gain
        % The Qlq and Rlq were defined in a previous part in this code
            %K = dlqr(Aa,Ba,Qlq,Rlq);
            %disp('LQR Gain K = '); disp(K);
        
        % Solving the DARE (Difference Algebraic Riccati Equation)
        P = 1e6*eye(n,n);
        for i = 1:1000
        P = Aa'*P*Aa -Aa'*P*Ba*inv(Ba'*P*Ba+Rlq)*Ba'*P*Aa+Qlq;
        traceP(i)=trace(P);
        end
        figure; plot(traceP); ylabel('Trace of P'); xlabel('Sample');
        K = ( Aa'*P*Ba*inv(Ba'*P*Ba+Rlq) )';
        disp('DARE-based LQR K ='); disp(K);
        
        
        
	%% Kalman Filter Gain
        % The Qkf and Rkf were defined in a previous part in this code
        %L = ( dlqr(Aa',Ca',Qkf,Rkf) )';
        %disp('Kalman Filter Gain L = '); disp(L);
        
        % Solving the Observer DARE (Difference Algebraic Riccati Equation)
        Po = 1e6*eye(n,n);
        for i = 1:1000
        Po = Aa*Po*Aa' -Aa*Po*Ca'*inv(Ca*Po*Ca'+Rkf)*Ca*Po*Aa'+Qkf;
        tracePo(i)=trace(Po);
        end
        figure; plot(tracePo); ylabel('Trace of Po'); xlabel('Sample');
        L = ( Aa*Po*Ca'*inv(Ca*Po*Ca'+Rkf) );
        disp('DARE-based KF L ='); disp(L);
    
        
%% Control system simulation
tfinal = 10; % total simulation time in seconds
N = round( tfinal/Ts ); % total number of samples
t = 0:Ts:N*Ts-Ts; % time vector for the plots

    r1(1:5)=0; r1(6:N)=1; % Cart position reference sequence (m)
    r2(1:5)=0; r2(6:N)=0; % Pendulum angle reference sequence (rad)

    % Disturbances at the output
      v = 1*wgn(1,N, 1e-3, 'linear'); % Noise
      v2(1:round(N/2))=0; v2(round(N/2)+1:N)=+0.2*1; % Load

for k = 1:1
    % Nominal model initial conditions
    x(:,:,k) = [0;0; deg2rad(2) ;0]; % Init angle at 90+2 degress
    y(:,:,k) = C*x(:,:,k);
    u(k) = 0;
    
    % Augmented model initial conditions (KF and LQR)
    xa(:,:,k) = [0;0;  0;0;0;0]; 
    ya(:,:,k) = Ca*xa(:,:,k);
    du(k) = 0;
    
        % Estimated state variables based on the augmented ones
        x1e(k)=0;
        x2e(k)=0;
        x3e(k)=0;
        x4e(k)=0;
end
for k = 2:N
    % State space simulation model
    x(:,:,k) = A*x(:,:,k-1) +B*u(k-1);
    y(:,:,k) = C*x(:,:,k) +[1;1]*v(k) +[1;0]*v2(k);
    
    
    % Augmented state space model
    xa(:,:,k) = Aa*xa(:,:,k-1) +Ba*du(k-1) +L*( y(:,:,k-1)-ya(:,:,k-1) ) ;
    ya(:,:,k) = Ca*xa(:,:,k);
        % Estimated state variables based on the augmented ones
            x1e(k)=x1e(k-1) +xa(3,k);
            x2e(k)=x2e(k-1) +xa(4,k);
            x3e(k)=x3e(k-1) +xa(5,k);
            x4e(k)=x4e(k-1) +xa(6,k);
    
    % State-feedback control law
    du(k) = K*( [r1(k); r2(k); 0;0;0;0] -xa(:,:,k) );
    
    % Passing du(k) through the discrete integrator
    u(k) = u(k-1) +du(k); %u(k)=r1(k);
end

figure; % I/O plots
subplot(311)
    plot(t,r1,'k'); hold;
    plot(t,y(1,:),'b'); plot(t,ya(1,:),':r');
    legend('r1(t)','y1(t)','y1a(t): estimated');
    ylabel('Cart Pos. (m)');
    
subplot(312)
    plot(t,r2,'k'); hold;
    plot(t,y(2,:),'b'); plot(t,ya(2,:),':r');
    legend('r2(t)','y2(t)','y2a(t): estimated');
    ylabel('Angle (rad)');

subplot(313)
    plot(t,u,'b');
    ylabel('Force (N)'); xlabel('Time (s)');
    
figure; % estimated state variables of the augmented estimator
subplot(321)
    plot(t,xa(1,:),'r');
    ylabel('xa_1(t)'); xlabel('Time (s)');
    title('State variables of the augmented estimator');

subplot(322)
    plot(t,xa(2,:),'r');
    ylabel('xa_2(t)'); xlabel('Time (s)');

subplot(323)
    plot(t,xa(3,:),'r');
    ylabel('xa_3(t)'); xlabel('Time (s)');

subplot(324)
    plot(t,xa(4,:),'r');
    ylabel('xa_4(t)'); xlabel('Time (s)');
    
subplot(325)
    plot(t,xa(5,:),'r');
    ylabel('xa_5(t)'); xlabel('Time (s)');

subplot(326)
    plot(t,xa(6,:),'r');
    ylabel('xa_6(t)'); xlabel('Time (s)');

    
figure; % estimated state variables based on the nominal system
subplot(411)
    plot(t,x(1,:),'b'); hold;
    plot(t,x1e,'r');
    ylabel('x_1(t): Cart Pos (m)');
    title('Estimated state variables based on the nominal system');

subplot(412)
    plot(t,x(2,:),'b'); hold;
    plot(t,x2e,'r');
    ylabel('x_2(t): Cart Vel (m/s)');
    
subplot(413)
    plot(t,x(3,:),'b'); hold;
    plot(t,x3e,'r');
    ylabel('x_3(t): Pend Ang (rad)');
    
subplot(414)
    plot(t,x(4,:),'b'); hold;
    plot(t,x4e,'r');
    ylabel('x_4(t): Pend Vel (red/s)');
    
    
%% Frequency response analysis in closed-loop

    % LQR analysis
    Tsen_lq = ss(Aa-Ba*K, Ba*K(1), Ca, Da, Ts);
    Ssen_lq = 1-Tsen_lq;
    figure; sigma(Tsen_lq);hold; sigma(Ssen_lq);
    title('LQR Sensitivity Functions');
    legend('|Tsen|','|Ssen|');
    
    
    % KF analysis
    Tsen_kf = ss(Aa-L*Ca, L, Ca, zeros(ny,ny), Ts);
    Ssen_kf = eye(ny,ny)-Tsen_kf;
    figure; sigma(Tsen_kf);hold; sigma(Ssen_kf);
    title('Kalmn Filter Sensitivity Functions');
    legend('|Tsen|','|Ssen|');
    
    % LQG analysis
    Adc = [Aa -Ba*K; L*Ca  Aa-Ba*K-L*Ca];
    Bdc = [Ba*K(1); Ba*K(1)];
    Cdc = [Ca zeros(ny,n)];
    Ddc = 0;
    
    Tsen_lqg = ss(Adc, Bdc, Cdc, Ddc, Ts);
    Ssen_lqg = [1;0]-Tsen_lqg;
    figure; sigma(Tsen_lqg);hold; sigma(Ssen_lqg);
    title('LQG Sensitivity Functions');
    legend('|Tsen|','|Ssen|');

    