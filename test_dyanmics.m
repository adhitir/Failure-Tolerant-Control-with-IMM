%% Testing the system dynamics and cable snapping

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

sim.dtdyn = 0.001;
sim.dt = 0.01;
sim.tf = 10;
sim.t_dyn = 0;

x0 = [0.4;0.43;0.0]; %m
x_dot0 = [0;0;0];
l0 = [0.1;0.1;0.1;0.1]; %m
theta0 = [0.3;0.3;0.3;0.3];

X_dyn(:,1) = [x0;x_dot0];

i=1;

model = m1;

%myVideo = VideoWriter('motionmodels3'); %open video file
%myVideo.FrameRate = 1000;  %can adjust this, 5 - 10 works well for me
%open(myVideo)

for t = 0:sim.dtdyn:sim.tf

    % Plot
    %model.plot_cables(X_dyn(1:3,i),l0)
    %%%%%%    
    
    if t == 2
        model = m2;
        sim.t_dyn = 0;
    end
    
    if t == 6
        model = m6;
        sim.t_dyn = 0;
    end
    sim.t_dyn = sim.t_dyn+sim.dtdyn;
       
    % Effect of joint inputs on the dynamic model.
    X_dyn_dot = model.fwd_dynamics(X_dyn(:,i),l0,theta0,sim.t_dyn);

    X_dyn(1:3,i+1) = X_dyn(1:3,i) + X_dyn_dot(1:3)*sim.dtdyn;
    X_dyn(4:6,i+1) = X_dyn_dot(1:3) + X_dyn_dot(4:6)*sim.dtdyn;

    k_vec(:,i) = model.k0current;

    %frame = getframe(gcf); %get frame
    %writeVideo(myVideo, frame);
    i = i+1;
end

%close(myVideo)

figure(1)
plot(0:sim.dtdyn:sim.tf,k_vec(:,1:i-1))
xlabel('time (s)');
ylabel('spring stiffness deterioration');
legend('k1','k2','k3','k4')

figure(2)
plot(0:sim.dtdyn:sim.tf,X_dyn(1:3,1:i-1))
xlabel('time (s)');
ylabel('states');
legend('x','y','\phi')
