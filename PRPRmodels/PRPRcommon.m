classdef PRPRcommon
   properties (Constant)
        base_length = 0.860; % mm
        base_span = 0.570; % mm

        phiA = 0;
        phiB = pi/2;
        phiC = pi;
        phiD = 3*pi/2;
        
        A1 = [PRPRcommon.base_length/2-PRPRcommon.base_span/2;0];
        A2 = [PRPRcommon.base_length/2+PRPRcommon.base_span/2;0];
        
        B1 = [PRPRcommon.base_length;PRPRcommon.base_length/2-PRPRcommon.base_span/2];
        B2 = [PRPRcommon.base_length;PRPRcommon.base_length/2+PRPRcommon.base_span/2];
        
        C1 = [PRPRcommon.base_length/2+PRPRcommon.base_span/2;PRPRcommon.base_length];
        C2 = [PRPRcommon.base_length/2-PRPRcommon.base_span/2;PRPRcommon.base_length];
        
        D1 = [0;PRPRcommon.base_length/2+PRPRcommon.base_span/2];
        D2 = [0;PRPRcommon.base_length/2-PRPRcommon.base_span/2];
        
        b = 0.16;
        h = 0.04;
        
        EA_O = [-PRPRcommon.b/2; -PRPRcommon.h/2];
        EB_O = [PRPRcommon.b/2; -PRPRcommon.h/2];
        EC_O = [PRPRcommon.b/2; PRPRcommon.h/2];
        ED_O = [-PRPRcommon.b/2; PRPRcommon.h/2];
        
        P = [0;0];
   end
   properties(Dependent)
       k0current
   end
   properties
       M
       D
       k0
       t_snap
   end

  methods 
      function obj = PRPRcommon(model)
           obj.M = [model.m*eye(2) [0;0]; [0 0] PRPRcommon.b*PRPRcommon.h*(PRPRcommon.b^2+PRPRcommon.h^2)/12];
           obj.D = [model.dx 0 0; 0 model.dy 0; 0 0 model.dz];
           obj.k0 = model.k0;
           obj.t_snap = model.t_snap;
      end

      function val = get.k0current(~)
           val = manage_k0;
      end

      function obj = set.k0current(obj,val)
          manage_k0(val);
      end
     
      function currentk0 = k0i(obj,t)
        % Let the cable fail over a period of 0.1 s
        currentk0 = obj.k0;
        if obj.fail > 0 && obj.fail <= 4
            currentk0(obj.fail) = obj.k0(obj.fail) - min(((obj.k0(obj.fail)-0)/obj.t_snap)*t,obj.k0(obj.fail));
        elseif obj.fail == 5
            currentk0(1) = min(obj.k0current(1),obj.k0(1) - min(((obj.k0(1)-0)/obj.t_snap)*t,obj.k0(1)));
            currentk0(3) = min(obj.k0current(3),obj.k0(3) - min(((obj.k0(3)-0)/obj.t_snap)*t,obj.k0(3)));
        elseif obj.fail == 6
            currentk0(2) = min(obj.k0current(2),obj.k0(2) - min(((obj.k0(2)-0)/obj.t_snap)*t,obj.k0(2)/2));
            currentk0(4) = min(obj.k0current(4),obj.k0(4) - min(((obj.k0(4)-0)/obj.t_snap)*t,obj.k0(4)/2));
        end
        obj.k0current = currentk0;
      end

      function X = fwd_kinematics(obj,X0,theta,l)
        
        %X0 = [0;0;0];
        Xmin = [0;0;-pi];
        Xmax = [1;1;pi];
        fnonlcon = @(X)ensure_positive_tension(X);
        
        X = fmincon(@(X)potential_energy(X), X0,[],[],[],[],Xmin,Xmax,fnonlcon);
        
        %k0 = k0i(obj,t)    
        function U = potential_energy(X)
            current_lp = obj.prismatic_length(X,l);
            k = obj.k0;
            switch obj.fail
                case {1,2,3,4}
                    k(obj.fail) = 0;
                case 5
                    k(1) = 0;
                    k(3) = 0;
                case 6
                    k(2) = 0;
                    k(4) = 0;
            end
            K = diag(k./theta);
            U = 0.5*(current_lp-theta)'*K*(current_lp-theta);
        end

        function [c,ceq] = ensure_positive_tension(X)
            current_lp = obj.prismatic_length(X,l);
            k = obj.k0;
            switch obj.fail
                case {1,2,3,4}
                    k(obj.fail) = 0;
                case 5
                    k(1) = 0;
                    k(3) = 0;
                case 6
                    k(2) = 0;
                    k(4) = 0;
            end
            K = diag(k./theta);
            c = -K*(current_lp-theta);
            ceq=[];
        end
      end
           
      function X_dyn_dot = fwd_dynamics(obj,X_dyn,l,theta,t)
        
        k0 = k0i(obj,t);    
        P_x = obj.structureMatrix1(X_dyn(1:3),l);
        l_p_x = obj.prismatic_length(X_dyn(1:3),l);
        K = diag(k0./l_p_x);
        x_dot = X_dyn(4:6);
        x_ddot = inv(obj.M)*(P_x*K*(l_p_x-theta)) - inv(obj.M)*obj.D*x_dot;
        
        X_dyn_dot = [x_dot;x_ddot];

      end
                
      function X_dyn_dot = fwd_dynamics_ekf(obj,X_dyn,l,theta,t)
        
        k0 = obj.k0;    
        if obj.fail > 0 && obj.fail <= 4
            k0(obj.fail) = 0;
        elseif obj.fail == 5
            k0(1) = 0;
            k0(3) = 0;
        elseif obj.fail == 6
            k0(2) = 0;
            k0(4) = 0; 
        end
        P_x = obj.structureMatrix1(X_dyn(1:3),l);
        l_p_x = obj.prismatic_length(X_dyn(1:3),l);
        K = diag(k0./l_p_x);
        x_dot = X_dyn(4:6);
        x_ddot = obj.M\(P_x*K*(l_p_x-theta)) - obj.M\obj.D*x_dot;
        
        X_dyn_dot = [x_dot;x_ddot];

      end
      
      function [A3,B3,C3,D3] = slider_positions(obj,l)
        A3 = obj.P + obj.A1 + [cos(obj.phiA) -sin(obj.phiA);sin(obj.phiA) cos(obj.phiA)]*[l(1);0];
        B3 = obj.P + obj.B1 + [cos(obj.phiB) -sin(obj.phiB);sin(obj.phiB) cos(obj.phiB)]*[l(2);0];
        C3 = obj.P + obj.C1 + [cos(obj.phiC) -sin(obj.phiC);sin(obj.phiC) cos(obj.phiC)]*[l(3);0];
        D3 = obj.P + obj.D1 + [cos(obj.phiD) -sin(obj.phiD);sin(obj.phiD) cos(obj.phiD)]*[l(4);0];  
      end

      function [EA_P,EB_P,EC_P,ED_P] = platform_coordinates(obj,X)
        pGo = obj.transformation_matrix(X);

        t1 = pGo*[obj.EA_O;1];
        t2 = pGo*[obj.EB_O;1];
        t3 = pGo*[obj.EC_O;1];
        t4 = pGo*[obj.ED_O;1];

        EA_P = t1(1:2);
        EB_P = t2(1:2);
        EC_P = t3(1:2);
        ED_P = t4(1:2);
      end

      function pGo = transformation_matrix(obj,X)
        phi_P = X(3);
        pGo = eye(3);
        pGo(1:2,1:2) = [cos(phi_P) -sin(phi_P); sin(phi_P) cos(phi_P)]; %rot_z
        pGo(1:2,3) = X(1:2);
      end

      function l_p = prismatic_length(obj,X,l)

        [EA_P,EB_P,EC_P,ED_P] = obj.platform_coordinates(X);       
        [A3,B3,C3,D3] = obj.slider_positions(l); 

        la = sqrt((A3(1)-EA_P(1))^2+(A3(2)-EA_P(2))^2);
        lb = sqrt((B3(1)-EB_P(1))^2+(B3(2)-EB_P(2))^2);
        lc = sqrt((C3(1)-EC_P(1))^2+(C3(2)-EC_P(2))^2);
        ld = sqrt((D3(1)-ED_P(1))^2+(D3(2)-ED_P(2))^2);

        l_p = [la,lb,lc,ld]';

      end
      
      function [theta1a,theta1b,theta1c,theta1d] = joint_angle(obj,X,l)
        [EA_P,EB_P,EC_P,ED_P] = obj.platform_coordinates(X);       
        [A3,B3,C3,D3] = obj.slider_positions(l); 

        va = [1,0]'; vb = [0,1]'; vc = [-1,0]'; vd = [0,-1]';
        ua = [(EA_P(1)-A3(1)),(EA_P(2)-A3(2))]'/(sqrt((EA_P(1)-A3(1))^2 + (EA_P(2)-A3(2))^2));
        ub = [(EB_P(1)-B3(1)),(EB_P(2)-B3(2))]'/(sqrt((EB_P(1)-B3(1))^2 + (EB_P(2)-B3(2))^2));
        uc = [(EC_P(1)-C3(1)),(EC_P(2)-C3(2))]'/(sqrt((EC_P(1)-C3(1))^2 + (EC_P(2)-C3(2))^2));
        ud = [(ED_P(1)-D3(1)),(ED_P(2)-D3(2))]'/(sqrt((ED_P(1)-D3(1))^2 + (ED_P(2)-D3(2))^2));

        theta1a = acos(max(min(dot(ua,va)/(norm(ua)*norm(va)),1),-1));
        theta1b = acos(max(min(dot(ub,vb)/(norm(ub)*norm(vb)),1),-1));
        theta1c = acos(max(min(dot(uc,vc)/(norm(uc)*norm(vc)),1),-1));
        theta1d = acos(max(min(dot(ud,vd)/(norm(ud)*norm(vd)),1),-1));
      end

      function Jw = structureMatrix(obj,X,l)
        phi_P = X(3);
        [theta1a,theta1b,theta1c,theta1d] = joint_angle(obj,X,l);

        gammaA = (phi_P - theta1a - obj.phiA);
        gammaB = (phi_P - theta1b - obj.phiB);
        gammaC = (phi_P - theta1c - obj.phiC);
        gammaD = (phi_P - theta1d - obj.phiD);

        % Pulling map:
        Jw_b = -[-cos(gammaA),                                        -cos(gammaB),                                        -cos(gammaC),                                       -cos(gammaD);
                 sin(gammaA),                                         sin(gammaB),                                         sin(gammaC),                                        sin(gammaD);
                 obj.EA_O(2)*cos(gammaA) + obj.EA_O(1)*sin(gammaA),   obj.EB_O(2)*cos(gammaB) + obj.EB_O(1)*sin(gammaB),   obj.EC_O(2)*cos(gammaC) + obj.EC_O(1)*sin(gammaC),  obj.ED_O(2)*cos(gammaD) + obj.ED_O(1)*sin(gammaD)];


        Jw = Jw_b;

      end
        
      function Jw = structureMatrix1(obj,X,l)
        [A3,B3,C3,D3] = obj.slider_positions(l); 
        
        phi_e = X(3);
        pRo = [cos(phi_e) -sin(phi_e); sin(phi_e) cos(phi_e)]; %rot_z
                
        la = A3(1:2) - X(1:2) - pRo*obj.EA_O;
        lb = B3(1:2) - X(1:2) - pRo*obj.EB_O;
        lc = C3(1:2) - X(1:2) - pRo*obj.EC_O;
        ld = D3(1:2) - X(1:2) - pRo*obj.ED_O;
        
        ua = [(la/norm(la));0];
        ub = [(lb/norm(lb));0];
        uc = [(lc/norm(lc));0];
        ud = [(ld/norm(ld));0]; 
        
        r1 = [pRo*obj.EA_O;0];
        r2 = [pRo*obj.EB_O;0];
        r3 = [pRo*obj.EC_O;0];
        r4 = [pRo*obj.ED_O;0];
        
        baxua = cross(r1,ua);
        bbxub = cross(r2,ub);
        bcxuc = cross(r3,uc);
        bdxud = cross(r4,ud);
                
        Jw = [ua(1:2), ub(1:2), uc(1:2), ud(1:2);
              baxua(3),bbxub(3),bcxuc(3),bdxud(3)];
        
      end

  end
end

function out=manage_k0(val) %Not a class method - can be class-related function
 persistent p
 if nargin
     p=val;
 end
 out=p;
end