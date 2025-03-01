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
options = odeset('RelTol',1e-6,'AbsTol',1e-6);

dp_ode = @(t, theta) [theta(3); theta(4);
    (-g * (2*m1 + m2) * sin(theta(1)) - m2 * g * sin(theta(1) - 2*theta(2)) - 2*sin(theta(1)-theta(2))*m2*(theta(4)^2 * L2 + theta(3)^2 * L1 * cos(theta(1)-theta(2)))) / (L1 * (2*m1 + m2 - m2 * cos(2*theta(1) - 2*theta(2))));
    (2 * sin(theta(1) - theta(2)) * (theta(3)^2 * L1 * (m1 + m2) + g * (m1 + m2) * cos(theta(1)) + theta(4)^2 * L2 * m2 * cos(theta(1)-theta(2)))) / (L2 * (2*m1 + m2 - m2 * cos(2*theta(1) - 2*theta(2))))];

tspan = [0:dt:runtime];
[t, sol] = ode45(dp_ode, tspan, [theta1_0, theta2_0, omega1_0, omega2_0], options);

x1 = L1 * sin(sol(:,1));
y1 = -L1 * cos(sol(:,1));
x2 = x1 + L2 * sin(sol(:,2));
y2 = y1 - L2 * cos(sol(:,2));

figure('Color', 'w');
subplot(2,1,1);
hold on;
h_line1 = plot([0, x1(1)], [0, y1(1)], 'k-', 'LineWidth', 2);
h_line2 = plot([x1(1), x2(1)], [y1(1), y2(1)], 'k-', 'LineWidth', 2);
h_dot1 = plot(x1(1), y1(1), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 10);
h_dot2 = plot(x2(1), y2(1), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 10);
title('Double Pendulum Motion');
axis equal; axis([-1.2*(L1+L2) 1.2*(L1+L2) -1.2*(L1+L2) 1.2*(L1+L2)]);
hold off;

subplot(2,1,2);
hold on;
h_theta1 = animatedline('Color', 'b', 'LineWidth', 2);
h_theta2 = animatedline('Color', 'r', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Angle (rad)');
title('Theta1 (Blue) and Theta2 (Red) over Time');
hold off;

for k = 1:length(t)

    set(h_line1, 'XData', [0, x1(k)], 'YData', [0, y1(k)]);
    set(h_line2, 'XData', [x1(k), x2(k)], 'YData', [y1(k), y2(k)]);
    set(h_dot1, 'XData', x1(k), 'YData', y1(k));
    set(h_dot2, 'XData', x2(k), 'YData', y2(k));

    addpoints(h_theta1, t(k), sol(k,1));
    addpoints(h_theta2, t(k), sol(k,2));

    drawnow;
end