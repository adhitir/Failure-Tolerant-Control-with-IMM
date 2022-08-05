
classdef twoPRPR < PRPRcommon
   properties
        fail
   end
   methods
      function obj = twoPRPR(model,fail)
           obj = obj@PRPRcommon(model);
           obj.fail = fail.fail;
       end
      function plot_cables(obj, X, l)      

        [EA_P,EB_P,EC_P,ED_P] = obj.platform_coordinates(X);
        [A3,B3,C3,D3] = obj.slider_positions(l); 
        
        figure(1);
        plot([obj.A1(1),obj.A2(1),obj.B1(1),obj.B2(1),obj.C1(1),obj.C2(1),obj.D1(1),obj.D2(1),obj.A1(1)],[obj.A1(2),obj.A2(2),obj.B1(2),obj.B2(2),obj.C1(2),obj.C2(2),obj.D1(2),obj.D2(2),obj.A1(2)],'b','LineWidth',1);
        hold on;
        
        % Plot ee
        plot([EA_P(1),EB_P(1),EC_P(1),ED_P(1),EA_P(1)],[EA_P(2),EB_P(2),EC_P(2),ED_P(2),EA_P(2)],'k','LineWidth',1);

        % Plot cables
        switch obj.fail 
            case 5
                plot([B3(1),EB_P(1)],[B3(2),EB_P(2)],'-.k','LineWidth',1.5);
                plot([D3(1),ED_P(1)],[D3(2),ED_P(2)],'-.k','LineWidth',1.5);
            case 6
                plot([A3(1),EA_P(1)],[A3(2),EA_P(2)],'-.k','LineWidth',1.5);
                plot([C3(1),EC_P(1)],[C3(2),EC_P(2)],'-.k','LineWidth',1.5);
        end
        axis equal
        hold off
      end
      
      function [theta,l,tau_pos] = joint_state(obj,X,X_dot,theta0,l0,dt,tau_pos_init)
          
        % New desired position of the end-effector
        x1_des = X + X_dot*dt;

        Jw = obj.structureMatrix1(X,l0);
        l_p_x = obj.prismatic_length(X,l0);
        [theta1a,theta1b,theta1c,theta1d] = joint_angle(obj,X,l0);
        B = [cos(theta1a) 0 0 0; 0 cos(theta1b) 0 0; 0 0 cos(theta1c) 0; 0 0 0 cos(theta1d)];
        
        switch obj.fail
            case 5
            l_p_x_dot = -Jw(1,:)'*X_dot(1);
            l_dot = -B\Jw(2,:)'*X_dot(2);
            case 6
            l_p_x_dot = -Jw(1,:)'*X_dot(1);
            l_dot = B\Jw(2,:)'*X_dot(2);
        end       
        
        l_dot(l_dot>0.05) = 0.05;
        l_dot(l_dot<-0.05) = -0.05;

        l = l0 + l_dot*dt;
        
        l(l>=obj.base_span) = obj.base_span;
        l(l<0) = 0;
        
        % Joint angles that will generate the desired pose and cable
        % tensions.

        l_p_x_new = l_p_x + l_p_x_dot*dt;
        k = obj.k0;
        switch obj.fail
            case 5
                tau_pos = [0;0.5;0;0.5];
            case 6
                tau_pos = [0.5;0;0.5;0];
        end

        theta = (k.*l_p_x_new)./(tau_pos + k);

        theta_dot = (theta - theta0)/dt;

        theta_dot(theta_dot>0.05) = 0.05;
        theta_dot(theta_dot<-0.05) = -0.05;

        theta = theta0 + theta_dot*dt; 
      
      end
      

      
      function [P_plus, X_hat_plus, likelihood] = EKF(obj,m,n,X_hat_plus,y_i,P_plus,thetai,li,t,dt,Q,R)

        % EKF
        f_xu = fwd_dynamics_ekf(obj,X_hat_plus,li,thetai,t);
        h_x = X_hat_plus(1:3);

        % Diff
        F = zeros(n,n);
        H = zeros(m,n);
        dx = f_xu*dt;
        for k = 1:n
            Xk_plus = X_hat_plus;
            Xk_plus(k) = Xk_plus(k) + dx(k);
            
            f_xuk = fwd_dynamics_ekf(obj,Xk_plus,li,thetai,t);
            h_xk = Xk_plus(1:3);

            Xk_minus = X_hat_plus;
            Xk_minus(k) = Xk_minus(k) - dx(k);

            f_xuk_ = fwd_dynamics_ekf(obj,Xk_minus,li,thetai,t);
            h_xk_ = Xk_minus(1:3);

            F(:,k) = (f_xuk - f_xuk_)/(2*dx(k));
            H(:,k) = (h_xk - h_xk_)/(2*dx(k));   

            F_discrete = eye(n) + F*dt;
        end

        % Propagation
        X_hat_minus = X_hat_plus + f_xu*dt;
        P_minus = F_discrete*P_plus*F_discrete' + Q;

        % Likelihood of measurement
        Ei = H*P_minus*H'+R;
        ei = y_i - h_x;
        likelihood = (1/sqrt(det(2*pi*Ei)))*exp(-0.5*ei'*inv(Ei)*ei);

        % Gain
        Kg = P_minus*H'*inv(H*P_minus*H'+R);

        % Update
        X_hat_plus = X_hat_minus + Kg*(y_i - h_x);
        P_plus = (eye(length(P_minus)) - Kg*H)*P_minus;
    end     
          
   end
end