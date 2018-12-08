clear
clc
g=9.8067;
k_t=0.3233;
k_b=0.4953;
R=5.2628;
L=0.11;
m_p=0.179;%no batteries
m_w=0.036;
r_w=0.016;
I_p=0.0014;
I_w=4.6*10^(-6);
A_2_1 = (g*L*m_p*(I_w + (m_p + m_w)*r_w^2))/(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2);
A_2_2 = -k_b*k_t*(I_w + r_w*(m_w*r_w + m_p*(L + r_w)))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
A_2_4 = -k_b*k_t*(I_w + r_w*(m_w*r_w + m_p*(L + r_w)))/(R*r_w*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
A_4_1 = (g*L^2*m_p^2*r_w^2)/(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2);
A_4_2 = -k_b*k_t*r_w*(I_p + L*m_p*(L + r_w))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
A_4_4 = -k_b*k_t*(I_p + L*m_p*(L + r_w))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
A=[0 1 0 0;
  A_2_1 A_2_2 0 A_2_4;
  0 0 0 1;
  A_4_1 A_4_2 0 A_4_4];
B_2 = -(k_t*(I_w + r_w*(m_w*r_w + m_p*(L + r_w))))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
B_4 = -(k_t*r_w*(I_p+ L*m_p*(L + r_w)))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
B = [0; B_2; 0; B_4];
C = eye(4);
D = zeros(4, 1);
%KLQR=[-500,-25,5,45];
KLQR=[-500,-25,5,47];%60 seconds
%KLQR=[-500,-25,10,150];
%KLQR=[-200,-25,10,150];
%KLQR=[7.7404, 47.163, -92.5425, -16.8075];
Q = diag([20000,20000,5,0.001]);
R = 100; 
%KLQR = lqr(A, B, Q, R);
encoder_counts=720;   % number of counts (if using quad encoding)
RPM_MAX = 170;        % spec sheet max RPM

% constants or conversion factors:
RADSEC2RPM = 60/(2*pi);       % radians/sec to RPM
RAD2C=encoder_counts/(2*pi);  % conversion from radians to counts
C2RAD=1/RAD2C; 
C2DEG=360/encoder_counts;
R2D=180/pi;
D2R=1/R2D;
g=9.8;

% hardware
Vsupply=9;  % 5 volts when on batter, 5 volts on USB
DCB2V=Vsupply/255;   % Duty cycle bits to PWM (volts)
V2DCB=1/DCB2V;       % Volts to duty cycle in bits.

% Poles placed  {-1065, -3.6+0.5, -0.4, -3.6-0.5i}.
%KLQR=[  7.7404, 47.163, -92.5425, -16.8075 ];   % 56mm wheel .381kg.
% Quite stable, not so much oscillation

% combine some constants:
ES=-C2DEG*D2R;

% automatic calibration of gyro offset
% discrete 1st order filter recurrence relation - 
% discrete-time implementation of a low-pass filter is 
% the exponentially-weighted moving average
% alpha=.5 (time contant is equal to sampling period)
% alpha<.5 (time constant is larger than sampling period) tau ~ TS/alpha

% .006 requires 800 samples to get to steardy state (see gyro filter design)
% 1st order system take 3Tau to get to .95 ss value, 3.9Tau, 98, 4.56tau 99
% tau ~ .0025/.006 = .4167 = 3*.4167 = 1.25 to get to .95  - we need .99
% tau ~ .0025/.006 = .4167 = .4167*.456 = 1.9sec to get to 99% ss
a_go = .006;  % alpha for initial gyro offset calibration (2 sec for ss)

% Start time for balancing - Gyro calibartion time:
%tstart=2;

%% MPU5060
%TS=.005 % Minimum simulation time
TS=0.005; % Simulation time with motor controller
%TS=.005 % fastest with default mpu5060 library settings... and filter set to DLPF set to 4
GyroS=250/32768;   % MPU5060 set to a maximum rate of 250 deg/second and 16-bit sampling
GS=GyroS*D2R;      % Convert from degrees to radians
tstart=.6;
return

% Low pass filter for output:
%a_ov = .7;  % on USB
%a_ov = .5; % on battery