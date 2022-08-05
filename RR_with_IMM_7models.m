%% IMM based FTC 
% Feedforward control

clc
clear
close all

addpath PRPRmodels

%% Common Model Parameters
model.m = 0.1; %kg
model.dx = 1;
model.dy = 1;
model.dz = 0.01;
model.k0 = [5;5;5;5]; %N
model.t_snap = 2; % Time taken for cable to fail

PRPRcommon(model);

%% Individual Model Parameters

model1.fail = 0;
model2.fail = 1;
model3.fail = 2;
model4.fail = 3;
model5.fail = 4;
model6.fail = 5;
model7.fail = 6;

m1 = fourPRPR(model,model1);
m2 = threePRPR(model,model2);
m3 = threePRPR(model,model3);
m4 = threePRPR(model,model4);
m5 = threePRPR(model,model5);
m6 = twoPRPR(model,model6);
m7 = twoPRPR(model,model7);

%% Initialize
sim.dt = 0.05;
sim.dt_dyn = 0.001;
sim.tf = 30cd;
sim.t = 0:sim.dt:sim.tf;
tn = length(sim.t);

x0 = [0.35;0.35;0.0]; %m
xdot0 = [0;0;0];

theta0 = [0.2;0.3;0.3;0.3];
l0 = [0;0;0;0]; %m
theta_dot = 0.001*ones(4,1);
l_dot = 0.001*ones(4,1);

%% Initialize the Controller
control.Kp = [0.8; 0.8; 0.0000];
control.Kd = [0.2; 0.2; 0.0000];
control.Kpl = 0.02;
control.A = 0.04; % Amplitude
control.periods = 1; % frequency
control.tf = sim.tf;
control.dt = sim.dt;
control.x0 = [0.35;0.35];
control.xf = [0.65;0.45];

c = controller(control);
%% Trajectory

[x1_traj,x2_traj] = c.sine_wave;

% hold on;
% plot(x1_traj(1,:),x1_traj(2,:),'--','Color','blue','LineWidth',1)
% m1.plot_cables(x1_traj(:,1),l0)

%% Initialize the IMM

model = m1;
% As model is 1, initial index is 1
ind = 1;

% Initialize the dynamic variables
n = 6; % number of states
m = 3; % number of measurements
X_dyn(:,1) = [x0;xdot0];
y(:,1) = [eye(3) zeros(3,3)]*X_dyn(:,1); % Measurements
upsilon = [0.0001;0.0001;0.0001;10;10;10];
q = 0.0001;
r = 0.00002;

et0 = zeros(3,1);
et01 = zeros(3,1);
et02 = zeros(3,1);
et03 = zeros(3,1);
et04 = zeros(3,1);
et05 = zeros(3,1);
et06 = zeros(3,1);
et07 = zeros(3,1);

tau_pos0 = [0.5;0.5;0.5;0.5];
tau_pos1 = [0.5;0.5;0.5;0.5];
tau_pos2 = [0.5;0.5;0.5;0.5];
tau_pos3 = [0.5;0.5;0.5;0.5];
tau_pos4 = [0.5;0.5;0.5;0.5];
tau_pos5 = [0.5;0.5;0.5;0.5];
tau_pos6 = [0.5;0.5;0.5;0.5];
tau_pos7 = [0.5;0.5;0.5;0.5];


x1_des = zeros(3,length(x1_traj)); %position
x2_des = zeros(3,length(x1_traj)); %velocity

x0 = model.fwd_kinematics(x0,theta0,l0);
X_dyn(:,1) = [x0;xdot0];
[theta(:,1),l(:,1),x1_kin(:,1),tau_pos0,et0] = c.kinematics_loop(m1,x0,theta0,l0,tau_pos0,et0,0,sim.dt);
 

% joint states:
theta(:,1) = theta0;
l(:,1) = l0;
theta1(:,1) = theta0;
l1(:,1) = l0;
theta2(:,1) = theta0;
l2(:,1) = l0;
theta3(:,1) = theta0;
l3(:,1) = l0;
theta4(:,1) = theta0;
l4(:,1) = l0;
theta5(:,1) = theta0;
l5(:,1) = l0;
theta6(:,1) = theta0;
l6(:,1) = l0;
theta7(:,1) = theta0;
l7(:,1) = l0;

x_kin1(:,1) = x0;
x_kin2(:,1) = x0;
x_kin3(:,1) = x0;
x_kin4(:,1) = x0;
x_kin5(:,1) = x0;
x_kin6(:,1) = x0;
x_kin7(:,1) = x0;

np = 7; % No of models;
weights_vec = zeros(np,length(x1_traj));
% weights(1,1:5/dt) = ones(1,5/dt);
% weights(2,5/dt+1:end) = ones(1,length(x1_traj)-5/dt);

R = eye(m)*r;
Q = upsilon*(0.00034^2)*upsilon';

P0 = eye(n)*0.1^2; %state covariance

X_hat = zeros(n,tn);
X_hat(:,1) = [0.430;0.40;0.0;0.0001;0.0001;0.0001];

X_hat_plus_bank = kron(ones(1,np),X_hat(:,1));
P_plus_bank(:,:,1) = P0;
P_plus_bank(:,:,2) = P0;
P_plus_bank(:,:,3) = P0;
P_plus_bank(:,:,4) = P0;
P_plus_bank(:,:,5) = P0;
P_plus_bank(:,:,6) = P0;
P_plus_bank(:,:,7) = P0;

mixed_initial_X = kron(ones(1,np),zeros(n,1));
mixed_initial_P = zeros(n,n,2);
                      
markov_transition_matrix=[0.70 0.05 0.05 0.05 0.05 0.05 0.05;
                          0.05 0.70 0.05 0.05 0.05 0.05 0.05;
                          0.05 0.05 0.70 0.05 0.05 0.05 0.05;
                          0.05 0.05 0.05 0.70 0.05 0.05 0.05;
                          0.05 0.05 0.05 0.05 0.70 0.05 0.05;
                          0.05 0.05 0.05 0.05 0.05 0.70 0.05;
                          0.05 0.05 0.05 0.05 0.05 0.05 0.70];
                      
likelihood = zeros(np,1);
weights = ones(np,1)/np;



%%
jd=1;

myVideo = VideoWriter('FTC'); %open video file
myVideo.FrameRate = 20;  %can adjust this, 5 - 10 works well for me
open(myVideo)

N = 7;
x_goal = [];
ti = 1;
sim.t_m = 0; % time from a mode change. 

for t = 0:sim.dt:sim.tf

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plant:
    % Cable snaps at the ~5s mark.
    if t == 10
        model = m2;
        t_m = 0;
    end
    
    if t == 15
        model = m6;
        t_m = 0;
    end
    
    for t_dyn = 0:sim.dt_dyn:sim.dt
        sim.t_m = sim.t_m+sim.dt_dyn;

        % Effect of joint inputs on the dynamic model.
        X_dyn_dot = model.fwd_dynamics(X_dyn(:,jd),l(:,ti),theta(:,ti),sim.t_m);

        X_dyn_dot = X_dyn_dot + upsilon*sqrt(q)*randn(1)*sqrt(sim.dt_dyn);

        X_dyn(1:3,jd+1) = X_dyn(1:3,jd) + X_dyn_dot(1:3)*sim.dt_dyn;
        X_dyn(4:6,jd+1) = X_dyn_dot(1:3) + X_dyn_dot(4:6)*sim.dt_dyn;

        %Generate measurements from plant model
        y(:,jd+1) = X_dyn(1:3,jd+1) + sqrt(r)*randn(m,1); 

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Measurement is fed to the IMM
        % IMM estimates the true state and weight vector

        %%% IMM %%%%
        % In the very first stage of the cycle, the weights and the estimates 
        % from the last cycle are mixed as per the Markov transition probabilities.
    
        % Mixed Weights
        for i = 1:np
            m_weights(:,i) = weights.*markov_transition_matrix(:,i);
        end

        normalizers = sum(m_weights);
    
        for i = 1:np
            for j = 1:np
                if normalizers(j) > 1e-20
                    m_weights(i,j) = m_weights(i,j) ./ normalizers(j);
                else
                    normalizers(j) = 0;
                    m_weights(i,j) = 0;
                end
            end
        end

        %keyboard

        % Mixed initial state
        for j = 1:np
            mixed_state = zeros(n, 1);
            for i = 1:np
                % mixed_state = mixed_State + w_ij * xi
                mixed_state = mixed_state + m_weights(i,j) * X_hat_plus_bank(:,i);
            end
            mixed_initial_X(:,j) = mixed_state;
        end
    
        %keyboard

        % mix initial state covariances
        for j = 1:np
            % for every filter, calculate the covaraince influence due
            % to every other filter

            % init zeros for covariance first
            mixed_cov = zeros(n);
            for i = 1:np
                % iterate through the wts of every filter

                %relative error from the mixed_init_mean
                err = X_hat_plus_bank(:,i) - mixed_initial_X(:,j);
                mixed_cov = mixed_cov +  m_weights(i,j)*(P_plus_bank(:,:,i) + err * err');
            end
            mixed_init_state_cov(:,:,j) = mixed_cov;
        end
        
        % EKF
        [P_plus_bank(:,:,1), X_hat_plus_bank(:,1), likelihood(1)] = m1.EKF(m,n,mixed_initial_X(:,1),y(:,jd+1),mixed_init_state_cov(:,:,1),theta(:,ti),l(:,ti),sim.t_m,sim.dt_dyn,Q,R); 
        [P_plus_bank(:,:,2), X_hat_plus_bank(:,2), likelihood(2)] = m2.EKF(m,n,mixed_initial_X(:,2),y(:,jd+1),mixed_init_state_cov(:,:,2),theta(:,ti),l(:,ti),sim.t_m,sim.dt_dyn,Q,R); 
        [P_plus_bank(:,:,3), X_hat_plus_bank(:,3), likelihood(3)] = m3.EKF(m,n,mixed_initial_X(:,3),y(:,jd+1),mixed_init_state_cov(:,:,3),theta(:,ti),l(:,ti),sim.t_m,sim.dt_dyn,Q,R); 
        [P_plus_bank(:,:,4), X_hat_plus_bank(:,4), likelihood(4)] = m4.EKF(m,n,mixed_initial_X(:,4),y(:,jd+1),mixed_init_state_cov(:,:,4),theta(:,ti),l(:,ti),sim.t_m,sim.dt_dyn,Q,R); 
        [P_plus_bank(:,:,5), X_hat_plus_bank(:,5), likelihood(5)] = m5.EKF(m,n,mixed_initial_X(:,5),y(:,jd+1),mixed_init_state_cov(:,:,5),theta(:,ti),l(:,ti),sim.t_m,sim.dt_dyn,Q,R); 
        [P_plus_bank(:,:,6), X_hat_plus_bank(:,6), likelihood(6)] = m6.EKF(m,n,mixed_initial_X(:,6),y(:,jd+1),mixed_init_state_cov(:,:,6),theta(:,ti),l(:,ti),sim.t_m,sim.dt_dyn,Q,R); 
        [P_plus_bank(:,:,7), X_hat_plus_bank(:,7), likelihood(7)] = m7.EKF(m,n,mixed_initial_X(:,7),y(:,jd+1),mixed_init_state_cov(:,:,7),theta(:,ti),l(:,ti),sim.t_m,sim.dt_dyn,Q,R); 

        if sum(likelihood) == 0 
            return
        end

        % Model Probability Update
        weights = weights.*likelihood;
        weights = weights/sum(weights); 

        % Conditional mean estimate
        P_plus = zeros(n);
        X_hat_plus(:,jd) = weights(1)*X_hat_plus_bank(:,1) + weights(2)*X_hat_plus_bank(:,2)  + weights(3)*X_hat_plus_bank(:,3) ...
            + weights(4)*X_hat_plus_bank(:,4)  + weights(5)*X_hat_plus_bank(:,5)  + weights(6)*X_hat_plus_bank(:,6)  + weights(7)*X_hat_plus_bank(:,7);
        for i = 1:np
            P_plus = P_plus + weights(i)*((X_hat_plus_bank(:,i) - X_hat_plus(:,jd))*(X_hat_plus_bank(:,i) - X_hat_plus(:,jd))' + P_plus_bank(:,:,i));
        end
        
        weights_vec(:,jd) = weights;
        P_cov(:,jd) =  diag(P_plus);
                
        %P_cov2(:,jd) = diag(chol(P_plus));
        
        jd = jd+1;

    end
       
    %keyboard
    % Weights determine which model is dominant.
        
    [theta1(:,ti+1),l1(:,ti+1),x_kin1(:,ti+1),tau_pos1(:,ti+1),et01] = c.kinematics_loop(m1,x1_traj(:,ti),theta(:,ti),l(:,ti),tau_pos1(:,ti),et01,sim.t_m,sim.dt);
    [theta2(:,ti+1),l2(:,ti+1),x_kin2(:,ti+1),tau_pos2(:,ti+1),et02] = c.kinematics_loop(m2,x1_traj(:,ti),theta(:,ti),l(:,ti),tau_pos2(:,ti),et02,sim.t_m,sim.dt);
    [theta3(:,ti+1),l3(:,ti+1),x_kin3(:,ti+1),tau_pos3(:,ti+1),et03] = c.kinematics_loop(m3,x1_traj(:,ti),theta(:,ti),l(:,ti),tau_pos3(:,ti),et03,sim.t_m,sim.dt);
    [theta4(:,ti+1),l4(:,ti+1),x_kin4(:,ti+1),tau_pos4(:,ti+1),et04] = c.kinematics_loop(m4,x1_traj(:,ti),theta(:,ti),l(:,ti),tau_pos4(:,ti),et04,sim.t_m,sim.dt);
    [theta5(:,ti+1),l5(:,ti+1),x_kin5(:,ti+1),tau_pos5(:,ti+1),et05] = c.kinematics_loop(m5,x1_traj(:,ti),theta(:,ti),l(:,ti),tau_pos5(:,ti),et05,sim.t_m,sim.dt);
    [theta6(:,ti+1),l6(:,ti+1),x_kin6(:,ti+1),tau_pos6(:,ti+1),et06] = c.kinematics_loop(m6,x1_traj(:,ti),theta(:,ti),l(:,ti),tau_pos6(:,ti),et06,sim.t_m,sim.dt);
    [theta7(:,ti+1),l7(:,ti+1),x_kin7(:,ti+1),tau_pos7(:,ti+1),et07] = c.kinematics_loop(m7,x1_traj(:,ti),theta(:,ti),l(:,ti),tau_pos7(:,ti),et07,sim.t_m,sim.dt);

    thetav = weights_vec(:,jd-1).*[theta1(:,ti+1)';theta2(:,ti+1)';theta3(:,ti+1)';theta4(:,ti+1)';theta5(:,ti+1)';theta6(:,ti+1)';theta7(:,ti+1)'];
    theta(:,ti+1) = sum(thetav,1);
    %theta(:,ti+1) = theta(:,ti) + 0.8*(sum(thetav,1)' -  theta(:,ti));
    
    lv = weights_vec(:,jd-1).*[l1(:,ti+1)';l2(:,ti+1)';l3(:,ti+1)';l4(:,ti+1)';l5(:,ti+1)';l6(:,ti+1)';l7(:,ti+1)'];
    l(:,ti+1) = sum(lv,1);
    %l(:,ti+1) = l(:,ti) + 0.1*(sum(lv,1)' -  l(:,ti));
    
    %keyboard;
    %%% Plotting
    %model.plot_cables(x_kin5(:,ti+1),l(:,ti))
    model.plot_cables(X_dyn(1:3,jd),l(:,ti))
    hold on;
    plot(x1_traj(1,1:ti),x1_traj(2,1:ti),'-','Color','blue','LineWidth',1);
    %plot(x1_des5(1,1:ti),x1_des5(2,1:ti),'-','Color','yellow','LineWidth',1);
    plot(X_dyn(1,1:jd),X_dyn(2,1:jd),'-','Color','red','LineWidth',1);
    hold off;
    xlabel('x (m)')
    ylabel('y (m)')
    %keyboard    

    %
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);

    %keyboard
    
    ti = ti+1;

end

%%
close(myVideo)


sig3_x=3*P_cov.^(0.5);

figure(2)
t_dyn = [0:1:jd-2]*0.00097;
plot(t_dyn,weights_vec,'LineWidth',1.5)
axis([0 sim.tf -0.2 1.2])
ylabel('Weights');
xlabel('Time (s)');
legend('No Failure','Cable A','Cable B','Cable C','Cable D','Cable A&C','Cable B&D');

figure(3)
t_dyn = [0:1:jd-2]*0.00097;
subplot(3,1,1)
plot(t_dyn,sig3_x(1,1:jd-1),'r',t_dyn,X_hat_plus(1,1:jd-1)-X_dyn(1,1:jd-1),'b',t_dyn,-sig3_x(1,1:jd-1),'r','LineWidth',1)
xlabel('Time (s)');
ylabel('Estimation Error (x)');
axis([0 15 -0.007 0.007 ])
subplot(3,1,2)
plot(t_dyn,sig3_x(2,1:jd-1),'r',t_dyn,X_hat_plus(2,1:jd-1)-X_dyn(2,1:jd-1),'b',t_dyn,-sig3_x(2,1:jd-1),'r','LineWidth',1)
xlabel('Time (s)');
ylabel('Estimation Error (y)');
axis([0 15 -0.007 0.007 ])
subplot(3,1,3)
plot(t_dyn,sig3_x(3,1:jd-1),'r',t_dyn,X_hat_plus(3,1:jd-1)-X_dyn(3,1:jd-1),'b',t_dyn,-sig3_x(3,1:jd-1),'r','LineWidth',1)
xlabel('Time (s)');
ylabel('Estimation Error (\phi)');
axis([0 15 -0.01 0.01 ])




figure(4)
t_dyn = [0:1:jd-1]*0.00097;
hold on;
plot(t_dyn,y(1,1:jd),'-r','LineWidth',1)
plot(t_dyn,X_dyn(1,1:jd),'-b','LineWidth',1)
%plot(t,x1_des1(1,1:length(x1_traj)),'-y','LineWidth',1)
plot(sim.t,x1_traj(1,1:length(x1_traj)),'-y','LineWidth',1)
xlabel('time (s)');
ylabel('States (x)');
title('True vs Measured State')
legend('x_m','x','x_{traj}')

figure(5)
subplot(3,1,1)
plot(sim.t,X_dyn(1,1:51:jd-1)-(x1_traj(1,1:length(x1_traj))),'-b','LineWidth',1)
xlabel('Time (s)');
ylabel('Trajectory Error(x)');
title('Trajectory Error correction during task recovery')
grid on;
subplot(3,1,2)
plot(sim.t,X_dyn(2,1:51:jd-1)-(x1_traj(2,1:length(x1_traj))),'-b','LineWidth',1)
xlabel('Time (s)');
ylabel('Trajectory Error(y)');
grid on;
subplot(3,1,3)
plot(sim.t,X_dyn(3,1:51:jd-1)-(x1_traj(3,1:length(x1_traj))),'-b','LineWidth',1)
xlabel('Time (s)');
ylabel('Trajectory Error(\phi)');
grid on;

figure(6)
d_error = sqrt((X_dyn(1,1:51:jd-1)-(x1_traj(1,1:length(x1_traj))))'.*(X_dyn(1,1:51:jd-1)-(x1_traj(1,1:length(x1_traj))))' + (X_dyn(2,1:51:jd-1)-(x1_traj(2,1:length(x1_traj))))'.*(X_dyn(2,1:51:jd-1)-(x1_traj(2,1:length(x1_traj))))');
plot(sim.t,d_error,'-b','LineWidth',1)
xlabel('Time (s)');
ylabel('Trajectory Error');
title('Trajectory Error correction during task recovery')
axis([0 20 -0.01 0.2])

figure(7)
subplot(4,1,1)
plot(sim.t,l(1,2:end),'-b','LineWidth',1.5)
xlabel('Time (s)');
ylabel('l_{a} (m)');
grid on;
subplot(4,1,2)
plot(sim.t,l(2,2:end),'-b','LineWidth',1.5)
xlabel('Time (s)');
ylabel('l_{b} (m)');
grid on;
subplot(4,1,3)
plot(sim.t,l(3,2:end),'-b','LineWidth',1.5)
xlabel('Time (s)');
ylabel('l_{c} (m)');
grid on;
subplot(4,1,4)
plot(sim.t,l(4,2:end),'-b','LineWidth',1.5)
xlabel('Time (s)');
ylabel('l_{d} (m)');
grid on;


figure(8)
plot(sim.t,l(1,2:end),'-','LineWidth',1.5)
hold on;
plot(sim.t,l(2,2:end),'-','LineWidth',1.5)
plot(sim.t,l(3,2:end),'-','LineWidth',1.5)
plot(sim.t,l(4,2:end),'-','LineWidth',1.5)
xlabel('Time (s)');
ylabel('Slider Position (m)');
legend('l_{a}', 'l_{b}', 'l_{c}', 'l_{d}');

