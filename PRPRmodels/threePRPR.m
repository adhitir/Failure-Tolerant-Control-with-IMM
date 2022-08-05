
classdef threePRPR < PRPRcommon
   properties
        fail
   end
   methods
       function obj = threePRPR(model,fail)
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
            case 1
                plot([B3(1),EB_P(1)],[B3(2),EB_P(2)],'-.k','LineWidth',1.5);
                plot([C3(1),EC_P(1)],[C3(2),EC_P(2)],'-.k','LineWidth',1.5);
                plot([D3(1),ED_P(1)],[D3(2),ED_P(2)],'-.k','LineWidth',1.5);
            case 2
                plot([A3(1),EA_P(1)],[A3(2),EA_P(2)],'-.k','LineWidth',1.5);
                plot([C3(1),EC_P(1)],[C3(2),EC_P(2)],'-.k','LineWidth',1.5);
                plot([D3(1),ED_P(1)],[D3(2),ED_P(2)],'-.k','LineWidth',1.5);
            case 3
                plot([A3(1),EA_P(1)],[A3(2),EA_P(2)],'-.k','LineWidth',1.5);
                plot([B3(1),EB_P(1)],[B3(2),EB_P(2)],'-.k','LineWidth',1.5);
                plot([D3(1),ED_P(1)],[D3(2),ED_P(2)],'-.k','LineWidth',1.5);
            case 4
                plot([A3(1),EA_P(1)],[A3(2),EA_P(2)],'-.k','LineWidth',1.5);
                plot([B3(1),EB_P(1)],[B3(2),EB_P(2)],'-.k','LineWidth',1.5);
                plot([C3(1),EC_P(1)],[C3(2),EC_P(2)],'-.k','LineWidth',1.5);
        end
        axis equal
        hold off
      end
       
      function [theta,l,tau_pos] = joint_state(obj,X,X_dot,theta0,l0,dt,tau_pos_init)

        x1_des = X + X_dot*dt;
          
        % Desired slider position for optimal manipulability ellipsoid
        l_des = obj.minimize_objective(x1_des,l0,theta0);

        % Realistic slider input that can be provided 
        delta_l = (l_des-l0)/dt;
        delta_l = max(min(0.005,delta_l),-0.005);
        l = l0 + delta_l;
                 
        Jw = obj.structureMatrix(x1_des,l);
        l_p_x = obj.prismatic_length(x1_des,l);
           

        switch obj.fail
            case 1
                Jw = Jw(1:2,[2,3,4]);
                tau_pos_init = tau_pos_init([2,3,4]);
            case 2
                Jw = Jw(1:2,[1,3,4]);
                tau_pos_init = tau_pos_init([1,3,4]);
            case 3
                Jw = Jw(1:2,[1,2,4]);
                tau_pos_init = tau_pos_init([1,2,4]);
            case 4
                Jw = Jw(1:2,[1,2,3]);
                tau_pos_init = tau_pos_init([1,2,3]);
        end
        
        % Determine corresponding joint angle
        f_o = [0;0];
        tau_min = [0.5;0.5;0.5];
        options = optimset('Display', 'off','LargeScale','on','Algorithm','interior-point');
        tau_pos = fmincon(@obj.minimize_tension, tau_pos_init ,[],[], Jw,-f_o, tau_min,[],[],options);

        switch obj.fail
            case 1
                tau_pos = [0;tau_pos];
            case 2
                tau_pos = [tau_pos(1);0;tau_pos(2:3)];
            case 3
                tau_pos = [tau_pos(1:2);0;tau_pos(3)];
            case 4
                tau_pos = [tau_pos;0];
        end
        
        k = obj.k0;
        %k(obj.fail) = 0;

        theta = (k.*l_p_x)./(tau_pos + k);

        theta_dot = (theta - theta0)/dt;

        theta_dot(theta_dot>0.05) = 0.05;
        theta_dot(theta_dot<-0.05) = -0.05;

        theta = theta0 + theta_dot*dt;
        
      end
      
      function tau = minimize_tension(obj,x)
          tau=norm(x,2);
      end
      
%       function X = fwd_kinematics(obj,X_init,theta,l)
%         
%         %X0 = [0;0;0];
%         Xmin = [0;0;-pi];
%         Xmax = [1;1;pi];
%         fnonlcon = @(X)ensure_positive_tension(X);
%         
%         [X,f] = fmincon(@(X)potential_energy(X), X_init,[],[],[],[],Xmin,Xmax,fnonlcon);
%         %k0 = k0i(obj,t)    
%         function U = potential_energy(X)
%             current_lp = obj.prismatic_length(X,l);
%             k = obj.k0;
%             k(obj.fail) = 0;
%             K = diag(k./theta);
%             U = 0.5*(current_lp-theta)'*K*(current_lp-theta);
%         end
% 
%        function [c,ceq] = ensure_positive_tension(X)
%             current_lp = obj.prismatic_length(X,l);
%             k = obj.k0;
%             k(obj.fail) = 0;
%             K = diag(k./theta);
%             c = -K*(current_lp-theta);
%             ceq=[];
%        end
%       end
            
      function l1 = minimize_objective(obj, x1i,l_init,theta)
        
        lb = [0;0;0;0];
        ub = [obj.base_span; obj.base_span; obj.base_span; obj.base_span];

        Aineq = [];
        bineq = [];
        Aeq = [];
        Beq = [];

        fnonlcon = @(l)ensure_wrench_closure(l);

        [l1,fval] = fmincon(@(l)objective_fn(l),l_init,Aineq,bineq,Aeq,Beq,lb,ub,fnonlcon);

        function objf = objective_fn(l)
            % Maximize wrench feasibility
            Jw = obj.structureMatrix(x1i,l);
            switch obj.fail
                case 1
                    Jw = Jw(:,[2,3,4]);
                case 2
                    Jw = Jw(:,[1,3,4]);
                case 3
                    Jw = Jw(:,[1,2,4]);
                case 4
                    Jw = Jw(:,[1,2,3]);
            end

            %A = Jw*Jw'; 
            %objf = -sqrt(det(A));  % MOM      
            z = null(Jw(1:2,:));
            if min(z) > 0
                objf = -min(z)/max(z);
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
            Ks = diag(k./theta);
            c = -Ks*(l_p_x-theta);
            
            % Ensure wrench closure
            Jw = obj.structureMatrix(x1i,l);
            switch obj.fail
                case 1
                    Jw = Jw(:,[2,3,4]);
                case 2
                    Jw = Jw(:,[1,3,4]);
                case 3
                    Jw = Jw(:,[1,2,4]);
                case 4
                    Jw = Jw(:,[1,2,3]);
            end
            z = null(Jw(1:2,:));
            z_bool = z>0;
            ceq = sum(z_bool) - 3;
        end
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