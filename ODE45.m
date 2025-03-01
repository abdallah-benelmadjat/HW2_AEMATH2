clc; clear; close all;

g = 9.81;   
m1 = 0.9;     
m2 = 0.2;   
L1 = 0.50;  
L2 = 0.30;  

theta1_0 = deg2rad(45);  
omega1_0 = deg2rad(0);   
theta2_0 = deg2rad(45);  
omega2_0 = deg2rad(0);   

runtime = 50;   
dt = 0.01;      
tspan = 0:dt:runtime;

options = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);

[t, Y] = ode45(@(t, y) double_pendulum_ODE(y, g, m1, m2, L1, L2), tspan, [theta1_0, theta2_0, omega1_0, omega2_0], options);

x1 = L1 * sin(Y(:, 1));
y1 = -L1 * cos(Y(:, 1));
x2 = x1 + L2 * sin(Y(:, 2));
y2 = y1 - L2 * cos(Y(:, 2));

figure('Color', 'w');
subplot(2,1,1);
hold on;
h_theta1 = animatedline('Color', 'b', 'LineWidth', 2);
h_theta2 = animatedline('Color', 'r', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Angle (rad)');
title('Theta1 (Blue) and Theta2 (Red) over Time');
hold off;

subplot(2,1,2);
hold on;
h_line1 = plot([0, x1(1)], [0, y1(1)], 'k-', 'LineWidth', 2);
h_line2 = plot([x1(1), x2(1)], [y1(1), y2(1)], 'k-', 'LineWidth', 2);
h_dot1 = plot(x1(1), y1(1), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
h_dot2 = plot(x2(1), y2(1), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
title('Double Pendulum Motion');
axis equal; axis([-1.2*(L1+L2) 1.2*(L1+L2) -1.2*(L1+L2) 1.2*(L1+L2)]);
hold off;

for k = 1:length(t)
    set(h_line1, 'XData', [0, x1(k)], 'YData', [0, y1(k)]);
    set(h_line2, 'XData', [x1(k), x2(k)], 'YData', [y1(k), y2(k)]);
    set(h_dot1, 'XData', x1(k), 'YData', y1(k));
    set(h_dot2, 'XData', x2(k), 'YData', y2(k));

    addpoints(h_theta1, t(k), Y(k,1));
    addpoints(h_theta2, t(k), Y(k,2));

    drawnow;
end

figure;
plot(t, rad2deg(Y(:, 1)), 'b', t, rad2deg(Y(:, 2)), 'r');
xlabel('Time (s)');
ylabel('Theta (degrees)');
legend('\theta_1', '\theta_2');
title('Double Pendulum Angles Over Time');
grid on;

function dydt = double_pendulum_ODE(y, g, m1, m2, L1, L2)
    theta1 = y(1); theta2 = y(2);
    omega1 = y(3); omega2 = y(4);

    delta = theta1 - theta2;
    den1 = (m1 + m2) * L1 - m2 * L1 * cos(delta)^2;
    den2 = (L2 / L1) * den1;

    domega1 = (-g*(2*m1 + m2)*sin(theta1) - m2*g*sin(theta1 - 2*theta2) - ...
               2*sin(delta)*m2*(omega2^2*L2 + omega1^2*L1*cos(delta))) / den1;

    domega2 = (2*sin(delta)*(omega1^2*L1*(m1+m2) + g*(m1+m2)*cos(theta1) + ...
              omega2^2*L2*m2*cos(delta))) / den2;

    dydt = [omega1; omega2; domega1; domega2];
end