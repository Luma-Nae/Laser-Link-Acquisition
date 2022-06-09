clc;
close all;

global Js K;

%% Simulation Setup
tspan = 0:0.1:10; % time vector spanning 10 hours



%% Satellite Parameters

Js1 = 1;
Js2 = 2;
Js3 = 5;
Js = diag([Js1, Js2, Js3]);

%% SC initial conditions


%      SC1
%ref frame initial orientation and angular velocity
q_SI = [ 1, 0, 0.1, 0]';
q_SI = q_SI/norm(q_SI);
omega_SI = [0; 0; 0]; %check this units



X0_1 = [q_SI; omega_SI]; %initial state vector


%      SC2

q_SI2 = [1, 0, 0.3, 0]';
%q_SI2 = q_SI2/norm(q_SI2);
omega_SI2 = [0; 0; 0];

q_r = quat2matrix(q_SI) * [q_SI2(1); -q_SI2(2:4)] ; % initial condition of relative quaternion state 
w_r = [0; 0; 0];


X0_r = [q_r; w_r]; %sc 2 initial state vector for the relative quaternion


%SC1 and SC2 initial state vector

X0_f = [X0_1 ; X0_r];





Xf = zeros(14, length(tspan)); %empty state vector spanned by time
Xf(:,1) = X0_f;

A = [zeros(3), 0.5*eye(3), zeros(3,6);
       zeros(3,12);
      zeros(3,9), 0.5*eye(3);
      zeros(3,12)];
B = [zeros(3), zeros(3); inv(Js), zeros(3); zeros(3,6); inv(Js), -inv(Js)];
Q = 10*eye(12);
R = 0.01*eye(6);
K = lqr(A,B,Q,R);


%% Simulation
figure(1);

for i = 2:length(tspan)
    opts = odeset('RelTol',1e-7,'AbsTol',1e-7);
     [~,y] = ode45(@SCAtDy, [tspan(i-1), tspan(i)], Xf(:, i-1), opts);
     Xf(:,i) = y(end,:)';  %output
     
     %normalizing quaternion
     Xf(1:4,i) = Xf(1:4,i) / norm(Xf(1:4,i));
     Xf(8:11,i) = Xf(8:11,i) / norm(Xf(8:11,i));
     
     R1 = q2R(Xf(1:4,i));
     T1 = RFrame(R1);

     Rr = q2R(Xf(8:11,i));
     Tr = RFrame(Rr);

     qrinv = [Xf(8,i);-Xf(9:11,i)];
     q2 = quat2matrix(qrinv)*Xf(1:4,i);
     R2 =  q2R(q2);
     T2 = RFrame(R2);

     clf;
     hold on;
     trplot(eye(4),'color','r'); %inertial reference frame plot
     trplot(T1);  %frame 1 plot
     trplot(T2, 'color', 'g'); %frame 2 plot
     trplot(Tr, 'color', 'm'); %relative quaternion plot
     grid on
     box on


     %drawing axes
     xlim([-1,1]);
     ylim([-1,1]);
     zlim([-1,1]);
     view(20,30);
     drawnow;
     
end


%% plots
figure; 
hold on
%SC 1 quaternion plot
plot (tspan, Xf(1,:)) 
plot (tspan, Xf(2,:)) 
plot (tspan, Xf(3,:)) 
plot (tspan, Xf(4,:)) 
title('SC1 attitude quaternion')

figure;
hold on
%SC2 quaternion plot
plot(tspan, Xf(8,:))
plot(tspan, Xf(9,:))
plot(tspan, Xf(10,:))
plot(tspan, Xf(11,:))

title('SC2 attitude quaternion')

figure; 
hold on
%SC1 angular velocity plot
plot (tspan, Xf(5,:))
plot (tspan, Xf(6,:)) 
plot (tspan, Xf(7,:)) 

title('SC1 angular velocity')

figure; 
hold on
%SC2 angular velocity plot
plot (tspan, Xf(12,:))
plot (tspan, Xf(13,:)) 
plot (tspan, Xf(14,:)) 

title('SC2 angular velocity')


%% Differential equation

function dX = SCAtDy(t, X)
    global Js K;
    
    %K1v2 = [K1(1:3,1:6),zeros(3,6)];

    q1 = X(1:4);
    w1 = X(5:7);
    qr = X(8:11);
    wr = X(12:14);

  %  input
   q1inv = [q1(1);-q1(2:4)];
   q_ref_s1 = [1.3, 0.2, 0.1, 0.1]';
   q_ref_s1 = q_ref_s1/norm(q_ref_s1);
   q_ss = quat2matrix(q1)*q_ref_s1;

   w_ss = w1;
   R21 = q2R(qr);
   R12 = inv(R21);

    tau = -K*[q_ss(2:4); w_ss; qr(2:4); wr];

    dX(1:4,1) = 1/2*quat2matrix(q1)*[0; w1];
    dX(5:7,1) = cross( -inv(Js) * w1,  w1) + inv(Js)*tau(1:3);

    dX(8:11,1) = 1/2*quat2matrix(qr)*[0; wr];
    dX(12:14,1) =cross(-inv(Js)*w_ss,w_ss) + inv(Js)*tau(1:3) ...
                        - cross(wr, (w_ss-wr)) ...
                        - cross(R21*inv(Js)*R12*(w_ss-wr), (w_ss-wr)) ...
                        -R21*inv(Js)*tau(4:6);
end


