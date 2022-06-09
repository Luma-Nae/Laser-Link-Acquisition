clc;
clear;
close all;

%% Parameters
Js1 = 1;
Js2 = 2;
Js3 = 5;
Js = diag([Js1, Js2, Js3]);

%% Defining symbolic variables

syms q1 q2 q3 q4 w1 w2 w3 qr1 qr2 qr3 qr4 wr1 wr2 wr3 t11 t12 t13 t21 t22 t23 w_ss

%states
w_ss = [w1; w2; w3]; 

w_r = [wr1; wr2; wr3]; 

%inputs
tau1 = [t11; t12; t13];

tau2 = [t21; t22; t23];

%matrices

S = [0, qr4, -qr3; -qr4, 0, qr2; qr3, -qr2, 0]; 
R12= eye(3)+ 2*S + 2*S^2;
R21= inv(R12); 


%operating points

w_ss_op = [0 , 0, 0]';


%%  Linearization
        %Spacecraft quaternion Linearization

A_qs = [ zeros(3), 0.5*eye(3), zeros(3,6)];

      %Spacecraft angular velocity linearization

f_w_ss = cross(inv(Js)* w_ss,  w_ss) - inv(Js)*tau1;

A_w_ss = jacobian(f_w_ss ,[ w_ss]); %linearized omega

%replacing operating points
A_w_ss = subs(A_w_ss, w_ss, w_ss_op);
A_w_ss = [zeros(3), A_w_ss, zeros(3,6)];


        %Relative quaternion Linearization
A_qr = [ zeros(3,9), 0.5*eye(3)];


        %Relative angular velocity Linearization

f_wr = cross(-inv(Js)*w_ss,w_ss) + inv(Js)*tau1 - cross(w_r, (-w_ss-w_r)) - cross(R21*inv(Js)*R12*(-w_ss-w_r), (-w_ss-w_r)) -R21*inv(Js)*tau2 ;
 
A_wr = jacobian(f_wr, [w_ss; qr2; qr3; qr4; w_r]);

A_wr = subs(A_wr, w_ss, [0; 0; 0]);
A_wr = subs(A_wr, qr2, [0]);
A_wr = subs(A_wr, qr3, [0]);
A_wr = subs(A_wr, qr4, [0]);
A_wr = subs(A_wr, w_r, [0; 0; 0]);
A_wr = subs(A_wr, tau1, [0;0;0]);
A_wr = subs(A_wr, tau2, [0;0;0]);

A_wr = [zeros(3), A_wr];
%connstructing a matrix

A = [A_qs ; A_w_ss; A_qr; A_wr];


        %B matrix

B_wr = jacobian(f_wr, [tau2]);
B_wr = subs(B_wr, tau2, [0; 0; 0]);
B_wr = subs(B_wr, qr2, [0]);
B_wr = subs(B_wr, qr3, [0]);
B_wr = subs(B_wr, qr4, [0]);


A = [zeros(3), 0.5*eye(3), zeros(3,6);
       zeros(3,12);
      zeros(3,9), 0.5*eye(3);
      zeros(3,12)];
B = [zeros(3), zeros(3); inv(Js), zeros(3); zeros(3,6); inv(Js), -inv(Js)];
Q = 10*eye(12);
R = 0.01*eye(6);
K = lqr(A,B,Q,R);

K1 = K(1:3,:);
K2 = K(4:6,:);








