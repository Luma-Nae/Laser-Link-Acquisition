%% Simulation setup

tspan = 0:0.1:365; 

%initial state vector
% x, y , z, x_dot, y_dot, z_dot
SC1_0 = [0.21052214 ;0.89889830; 0.40786493; -0.016670363; 0.0032193213;  0.0013952001];

SC2_0 = [0.19085381; 0.89912108; 0.38098013; -0.016938360; 0.0030315238; 0.0015874632];

SC3_0 = [0.22353943; 0.89276563; 0.37821601; -0.016880581; 0.0034893057; 0.0012394541];


X = zeros(18, length(tspan)); %empty state vector
X(:,1) = [SC1_0; SC2_0; SC3_0];

%% SIMULATION
figure;
hold on
for i = 2:length(tspan)
    opts = odeset('RelTol',1e-7,'AbsTol',1e-7);
    [~,y] = ode45(@OrbitD, [tspan(i-1), tspan(i)], X(:, i-1), opts);
    X(:,i) = y(end,:);
    if mod(i,300)==0
        scatter3(X(1,i),X(2,i),X(3,i),40,'filled','r')
        scatter3(X(7,i),X(8,i),X(9,i),40,'filled','g')
        scatter3(X(13,i),X(14,i),X(15,i),40,'filled', 'b')
    end
end


%% PLOTS

plot3 (X(1,:),X(2,:),X(3,:),'r')
plot3 (X(7,:),X(8,:),X(9,:),'g')
plot3 (X(13,:),X(14,:),X(15,:),'b')
axis equal
grid on
box on

%xlim([X(1,1)-0.1 X(1,1)+0.1]);
%ylim([X(2,1)-0.1 X(2,1)+0.1]);
%zlim([X(3,1)-0.1 X(3,1)+0.1]);

xlabel('X [AU]');
ylabel('Y [AU]');
zlabel('Z [AU]')
title('Initial Orbital Position')

%% Function
function dX = OrbitD(t, X)

    G = 6.67428*10^(-11) * 86400^2/(1.4957*10^(11) )^3; %gravitational constant in AU^3/kg*d^2
    Msun = 1988500*10^(24);

    r1 = X(1:3);
    v1 = X(4:6);

    r2 = X(7:9);
    v2 = X(10:12);

    r3 = X(13:15);
    v3 = X(16:18);

    dX(1:3,1) = v1;
    dX(4:6,1) = -(G*Msun)/norm(r1)^3*r1; 

    dX(7:9,1) = v2;
    dX(10:12,1) = -(G*Msun)/norm(r2)^3*r2; 

    dX(13:15,1) = v3;
    dX(16:18,1) = -(G*Msun)/norm(r3)^3*r3; 
end



