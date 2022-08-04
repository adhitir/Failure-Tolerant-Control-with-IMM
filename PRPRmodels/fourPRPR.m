
classdef fourPRPR < PRPRcommon
   properties
        fail 
   end
   methods
       function obj = fourPRPR(model,fail)
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
        plot([A3(1),EA_P(1)],[A3(2),EA_P(2)],'-.k','LineWidth',1.5);
        plot([B3(1),EB_P(1)],[B3(2),EB_P(2)],'-.k','LineWidth',1.5);
        plot([C3(1),EC_P(1)],[C3(2),EC_P(2)],'-.k','LineWidth',1.5);
        plot([D3(1),ED_P(1)],[D3(2),ED_P(2)],'-.k','LineWidth',1.5);

        
        axis equal
        hold off
      end
                        
      function [theta,l,tau_pos] = joint_state(obj,X,X_dot,theta0,l0,dt,tau_pos_init)

         
        % New desired position of the end-effector
        X_new = X + X_dot*dt;
        
        % Convert spatial velocity to body velocity
        phi_e = X(3);
        pRo = [cos(phi_e) -sin(phi_e); sin(phi_e) cos(phi_e)]; %rot_z
        p = X(1:2);      
        Adgpo = [pRo [p(2); -p(1)]; 0 0 1]; 
        X_dot_b = Adgpo\X_dot;
        
        
        % Desired slider position for optimal manipulability ellipsoid
        l_des = obj.minimize_objective(X_new,l0,theta0);
        % Realistic slider input that can be provided 
        l_dot = (l_des-l0)/dt;
        
        l_dot(l_dot>0.05) = 0.05;
        l_dot(l_dot<-0.05) = -0.05;

        l = l0 + l_dot*dt;
          
        % New desired prismatic length
        Jw = obj.structureMatrix1(X,l);
        l_p_x = obj.prismatic_length(X,l);
        
        l_p_x_dot = -Jw'*X_dot;
        l_p_x_new = l_p_x + l_p_x_dot*dt;
        
        Jw_new = obj.structureMatrix1(X_new,l);
        %l_p_x_new = obj.prismatic_length(X_new,l);
        
        % Is wrench closure?
%         z = null(Jw_new);
%         z_bool = z>0;
%         ceq = sum(z_bool) - 4;
%         if ceq ~= 0
%             keyboard
%         end
                
        % Desired tension in the new position
        f_o = [0;0;0];%0.1*[X_dot_b(1:2);0]; % desired wrench on ee is applied in the direction of the desired velocity.
        
        tau_min = [0.5;0.5;0.5;0.5];
        options = optimset('Display', 'off','LargeScale','on','Algorithm','interior-point');
        tau_pos = fmincon(@obj.minimize_tension, tau_pos_init ,[],[], Jw_new,-f_o, tau_min,[],[],options);
               
        k = obj.k0;

        theta = (k.*l_p_x_new)./(tau_pos + k);

        theta_dot = (theta - theta0)/dt;

        theta_dot(theta_dot>0.05) = 0.05;
        theta_dot(theta_dot<-0.05) = -0.05;

        theta = theta0 + theta_dot*dt;

%         Ks = k./l_p_x_new;
%         
%         % Joint angles that will generate the desired pose and cable
%         % tensions.
%         theta_des = l_p_x_new - tau_pos./Ks;
%         
%         theta =  theta_des;% - 0.5*(theta_des - theta0);
%         theta(isinf(theta)) = -100;
%         theta(isnan(theta)) = -100;
        
        %theta = l_p_x_new;
      end
    
      function tau = minimize_tension(obj,x)
          tau=norm(x,2);
      end
      
%       function X = fwd_kinematics(obj,X0,theta,l)
%         
%         %X0 = [0;0;0];
%         Xmin = [0;0;-pi];
%         Xmax = [1;1;pi];
%         fnonlcon = @(X)ensure_positive_tension(X);
%         
%         [X,f] = fmincon(@(X)potential_energy(X), X0,[],[],[],[],Xmin,Xmax,fnonlcon);
%         %k0 = k0i(obj,t)    
%         function U = potential_energy(X)
%             current_lp = obj.prismatic_length(X,l);
%             k = obj.k0;
%             %K = diag(k./current_lp);
%             K = diag(k./theta);
% 
%             U = 0.5*(current_lp-theta)'*K*(current_lp-theta);
%         end
% 
%        function [c,ceq] = ensure_positive_tension(X)
%             current_lp = obj.prismatic_length(X,l);
%             k = obj.k0;
%             %K = diag(k./current_lp);
%             K = diag(k./theta);
%             c = -K*(current_lp-theta);
%             
%             % Ensure wrench closure
%             Jw = obj.structureMatrix(X,l);
%             z = null(Jw);
%             z_bool = z>0;
%             ceq = sum(z_bool) - 4;
%        end
%       end
      
      
      function l1 = minimize_objective(obj, x1i,l_init,theta)
      % Maximize the manipulability ellipsoid.
        lb = [0;0;0;0];
        ub = [obj.base_span; obj.base_span; obj.base_span; obj.base_span];

        Aineq = [];
        bineq = [];
        Aeq = [];
        Beq = [];

        fnonlcon = @(l)ensure_wrench_closure(l);

        [l1,fval] = fmincon(@(l)objective_fn(l),l_init,Aineq,bineq,Aeq,Beq,lb,ub,fnonlcon);

        function objf = objective_fn(l)
            Jw = obj.structureMatrix(x1i,l);
%             A = Jw*Jw'; 
%             objf = -sqrt(det(A));  % MOM
            z = null(Jw);
            if min(z) > 0
                objf = -min(z)/max(z); %Sensitivity
            elseif max(z) < 0
                objf = -max(z)/min(z);
            else
                objf = 0;
            end
        end

        function [c,ceq] = ensure_wrench_closure(l)
            % Ensure positive tensions
            l_p_x = obj.prismatic_length(x1i,l);
            k = obj.k0;
            %Ks = diag(obj.k0./l_p_x);
            Ks = diag(k./theta);
            c = -Ks*(l_p_x-theta);
            
            % Ensure wrench closure
            Jw = obj.structureMatrix(x1i,l);
            z = null(Jw);
            z_bool = z>0;
            ceq = sum(z_bool) - 4;
        end
      end
             
      function [P_plus, X_hat_plus, likelihood] = EKF(obj,m,n,X_hat_plus,y_i,P_plus,thetai,li,t,dt,Q,R)

        % EKF
        g_xu = fwd_dynamics_ekf(obj,X_hat_plus,li,thetai,t);
        h_x = X_hat_plus(1:3);

        % Propagation
        % Discretize
        X_hat_minus = X_hat_plus + g_xu*dt;

        % To linearize, find dg(X)/d(X)
        G = zeros(n,n);
        H = zeros(m,n);
        dx = g_xu*dt;
        for k = 1:n
            Xk_plus = X_hat_plus;
            Xk_plus(k) = Xk_plus(k) + dx(k);
            
            g_xuk = fwd_dynamics_ekf(obj,Xk_plus,li,thetai,t);
            h_xk = Xk_plus(1:3);

            Xk_minus = X_hat_plus;
            Xk_minus(k) = Xk_minus(k) - dx(k);

            g_xuk_ = fwd_dynamics_ekf(obj,Xk_minus,li,thetai,t);
            h_xk_ = Xk_minus(1:3);

            G(:,k) = (g_xuk - g_xuk_)/(2*dx(k));
            H(:,k) = (h_xk - h_xk_)/(2*dx(k));   
        end

        F_discrete = eye(n) + G*dt;

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