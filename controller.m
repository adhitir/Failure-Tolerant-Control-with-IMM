classdef controller
    properties 
        Kp
        Kd
        Kpl
        A % Amplitude
        periods % frequency
        tf
        dt
        x0
        xf
    end
    methods
        function obj = controller(control)
           obj.Kp = control.Kp;
           obj.Kd = control.Kd;
           obj.Kpl = control.Kpl;
           obj.A = control.A;
           obj.periods = control.periods;
           obj.tf = control.tf;
           obj.dt = control.dt;
           obj.x0 = control.x0;
           obj.xf = control.xf;
       end

        function [x1,x2] = sine_wave(obj)
            
            x = obj.xf-obj.x0;
            dist = norm(obj.xf-obj.x0);
            t = 0:obj.dt:obj.tf;

            theta = atan2(x(2),x(1));
            R_theta = [cos(theta) -sin(theta); sin(theta) cos(theta)];

            omega = obj.periods*2*pi;
            
            X_des = obj.x0 + R_theta*[(dist/obj.tf)*t; obj.A*sin(omega*t/obj.tf)];

            X_dot_des = R_theta*[(dist/obj.tf)*ones(1,length(t)); obj.A*(omega/obj.tf)*cos((omega/obj.tf)*t)];

            phi_traj = atan2(X_dot_des(2,:),X_dot_des(1,:));

            X_ddot_des = R_theta*[zeros(1,length(t)); -obj.A*(omega/obj.tf)*(omega/obj.tf)*sin((omega/obj.tf)*t)];

            phi_dot_traj = (X_dot_des(1,:).*X_ddot_des(2,:) - X_dot_des(2,:).*X_ddot_des(1,:))./(X_dot_des(1,:).^2);

            x1 = [X_des;phi_traj];
            x2 = [X_dot_des;phi_dot_traj];
        end

        function [theta,l,x1,tau_pos,et0] = kinematics_loop(obj,model,x1_traji,theta,l,tau_pos_init,et0,t,dt)

            x1 = model.fwd_kinematics(x1_traji,theta,l);
    
            et = (x1_traji-x1);   
            x1_des = x1 + obj.Kp.*et + obj.Kd.*(et-et0);
            x2_des = (x1_des - x1)/dt;
            et0 = et;
            
            
            % theta for tracking trajectory
            [theta,l,tau_pos] = model.joint_state(x1,x2_des,theta,l,dt,tau_pos_init);
        end
        
    end
end
        
