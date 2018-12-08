close all
clear
clc
%% Parameters
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
%% Step 1&2 
A_2_1 = (g*L*m_p*(I_w + (m_p + m_w)*r_w^2))/(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2);
A_2_2 = -k_b*k_t*(I_w + r_w*(m_w*r_w + m_p*(L + r_w)))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
A_2_4 = -k_b*k_t*(I_w + r_w*(m_w*r_w + m_p*(L + r_w)))/(R*r_w*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
A_4_1 = (g*L^2*m_p^2*r_w^2)/(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2);
A_4_2 = -k_b*k_t*r_w*(I_p + L*m_p*(L + r_w))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
A_4_4 = -k_b*k_t*(I_p + L*m_p*(L + r_w))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
A=[0 1 0 0;
  A_2_1 A_2_2 0 A_2_4;
  0 0 0 1;
  A_4_1 A_4_2 0 A_4_4]
B_2 = -(k_t*(I_w + r_w*(m_w*r_w + m_p*(L + r_w))))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
B_4 = -(k_t*r_w*(I_p+ L*m_p*(L + r_w)))/(R*(I_w*(I_p + L^2*m_p) + (L^2*m_p*m_w + I_p*(m_p + m_w))*r_w^2));
B = [0; B_2; 0; B_4]
C = eye(4)
D = zeros(4, 1)
%% Step 3
sys = ss(A,B,C,D);
[num, den] = ss2tf(A, B, C, D);

%% Step 4
char_polyA=poly(A)
eigA=eig(A)
sys_cl = ss(A, B, C, D);
figure
pzmap(sys_cl);
%% Step 5
if(eigA < 0)
 fprintf('The system is asymptotically stable.\n')
else
 fprintf('The system is not asymptotically stable.\n')
end

if(eigA <= 0)
 fprintf('The system is marginally stable.\n')
else
 fprintf('The system is not marginally stable.\n')
end

%% Step 6
poles = pole(sys)
if(poles < 0 )
 fprintf('The system is BIBO stable.\n')
else
 fprintf('The system is not BIBO stable.\n')
end
%% Step 7
ctrb_matrix = ctrb(A,B)
ctrb_rank = rank(ctrb_matrix)
if(ctrb_rank<length(A))
 fprintf('The system is not controllable.\n')
else
 fprintf('The system is controllable.\n')
end
%% Step 8
obsv_matrix = obsv(A,C)
obsv_rank = rank(obsv_matrix)
% check the rank of the observability matrix
if(obsv_rank < length(A))
 fprintf('The system is not observable.\n\n')
else
 fprintf('The system is observable.\n\n')
end
%% Step 9
[num,den] = ss2tf(A,B,C,D);
AC=[-den(2),-den(3),-den(4),-den(5);  
    1 0 0 0;
    0 1 0 0;
    0 0 1 0];
BC=[1;0;0;0];
CC=[num(:,2) num(:,3) num(:,4) num(:,5)]; %???
DC=[0;0;0;0];
AO=[0 0 0 -den(5);
    1 0 0 -den(4);
    0 1 0 -den(3);
    0 0 1 -den(2)];
BO=[num(:,5)';  %???
    num(:,4)';
    num(:,3)';
    num(:,2)'];
CO=[0 0 0 1];
DO=0;
CCF = ss(AC,BC,CC,DC)
OCF = ss(AO,BO,CO,DO) 
%% Step 10&11
xini=[0; 0; 0; 0];
xhatini = xini;
tspan = 0:0.1:10;
pole1 = [-1,-2,-3+1i,-3-1i]; % target pole
L = place(transpose(A), transpose(C), pole1);
L = transpose(L);
sim('state_estimator',tspan(end))
% figure;
% T1 = t(1:209);
% T2 = t(2:210);
% delta_T = T2-T1
% plot(delta_T)
figure(1);
% title('State Estimator');
subplot(3,1,1); 
plot(t, x, 'LineWidth', 2.5); 
grid on;
xlabel('time (sec)'); 
legend( '\alpha', 'd\alpha', 'x','dx');
title('Closed-loop state estimator for the open-loop system')
% hold on;   
subplot(3,1,2); 
plot(t, x_hat, '-*'); 
grid on;   
xlabel('time (sec)'); 
legend( '\alpha_o_b_s', 'd\alpha_o_b_s', ...    
    'x_o_b_s','dx_o_b_s');
subplot(3,1,3); 
plot(t, y, '-.'); 
grid on;   
xlabel('time (sec)'); 
legend('y_1','y_2','y_3','y_4'); % Second pole location

%% Step 12-14
% close loop
pole_desired = [-1,-2,-3+1i,-3-1i]; % target pole 
K = place(A, B, pole_desired)
ACL = A - B*K; % Matrix A in the closed-loop system
char_polyACL=poly(ACL)
eigACL = eig(ACL)
sys_cl = ss(ACL, B, C, D);
figure
pzmap(sys_cl);
% If all the eigenvalues are located in LHP, the system is stable.
if(eigACL < 0)
 fprintf('The system is asymptotically stable.\n')
else
 fprintf('The system is not asymptotically stable.\n')
end
sim('close_loop',tspan(end));
figure(2);
% subplot(2,1,1); 
% plot(t, x, 'LineWidth', 2.5); 
title('Feedback Control');
% hold on;   
% xlabel('time (sec)'); 
% legend('a', 'da', 'x','dx');
% subplot(2,1,2); 
plot(t, y, '-.',  'LineWidth', 1.5); 
grid on;   
xlabel('time (sec)'); 
legend('\alpha', 'd\alpha', 'x','dx'); % Second pole location
title('Feedback Control');
%% Step 15-16
sim('combine',tspan(end));
figure(3);
title('Feedback Control using State Estimator');
subplot(2,1,1); 
plot(t, x, 'LineWidth', 1.5); 
grid on;   
xlabel('time (sec)'); 
legend('\alpha', 'd\alpha', 'x','dx');
subplot(2,1,2); 
plot(t, x_hat, '-*'); 
grid on;   
xlabel('time (sec)'); 
legend('\alpha_o_b_s', 'd\alpha_o_b_s', ...    
    'x_o_b_s','dx_o_b_s');
figure(4);
plot(t, x-x_hat, 'LineWidth', 1.5); 
title('State Estimator Error')
legend('\alpha', 'd\alpha', 'x','dx');