%% Initial Conditions
clear; close all; clc;
x0 = 5000; v0 = -100; m = 1000; Tmax = 30250; g = 9.81;

%% (a) Case 1: Max Thrust --> Free Fall
% Case 1
a = 1/2 * g * (1-m*g/Tmax);
b = m * v0 * g / Tmax;
c = x0 - 1/2 * m * v0^2 / Tmax;

syms x
eqn = a * x^2 + b * x + c == 0;
x_sol = solve(eqn);
t_f = double(vpa(x_sol,4));
t_trans = m*g/Tmax * t_f - m*v0/Tmax;

disp(['Final Time: ',num2str(t_f')]);
disp(['Thrust Conversion Time: ',num2str(t_trans')])
%% (a) Case 2: Free Fall --> Max Thrust
a = 1/2 * g * (m*g/Tmax - 1);
b = (1 - m*g/Tmax) * v0;
c = x0 + 1/2 * m * v0^2 / Tmax;

syms x
eqn = a * x^2 + b * x + c == 0;
x_sol = solve(eqn);
t_f = double(vpa(x_sol,4));
t_trans = m*v0/Tmax + t_f * (1 - m*g/Tmax);

disp(['Final Time: ',num2str(t_f')]);
disp(['Thrust Conversion Time: ',num2str(t_trans')])


%% Simulation with RK4
h = 0.001;

% Assumed scenario 2 
% Free Fall --> Thrust
for i=1:2
    if isreal(t_f(i)) && t_f(i) > 0 && t_trans(i) > 0
        X = [x0; v0];
        T = 0;
        t = 0;
        while t < t_f(i)
            X_curr = X(:,end);
            k1 = func(X_curr,t,t_trans(i),m,Tmax);
            k2 = func(X_curr+h*k1/2,t+h/2,t_trans(i),m,Tmax);
            k3 = func(X_curr+h*k2/2,t+h/2,t_trans(i),m,Tmax);
            k4 = func(X_curr+h*k3,t+h,t_trans(i),m,Tmax);
            X_next = X_curr + 1/6 * (k1 + 2*k2 + 2*k3 + k4) * h;
            X = [X, X_next];
            t = t + h;
            T = [T, t];
        end
    end
end

%% Figures
figure(1); axis on; grid on;
plot(T,X(1,:),'r.');
title('Rocket Position');
xlabel('Time(s)'); ylabel('Altitude(m)');

figure(2); axis on; grid on;
plot(T,X(2,:),'r.');
title('Rocket Velocity');
xlabel('Time(s)'); ylabel('Velocity(m/s)');

figure(3); axis on; grid on;
plot(X(1,:),X(2,:),'k.')
title('Phase Diagram')
xlabel('Position(X)'); ylabel('Velocity(V)');

height = 300;
width = 10;

figure(4); xlabel('X'); ylabel('Altitude(m)');
title('Rocket Landing Simulation');
axis([-50 50 0 x0])
r = rectangle('EdgeColor','k','FaceColor','b');
flag = false;
p = text(25,1/2*x0,'');
p.FontSize = 14;
% pause(10)

for i=1:50:size(X,2)
    % If free fall: blue
    % If thrust: red
    r.Position = [-1/2*width X(1,i) width height];
    if ~flag && T(i) > t_trans(2)
        r.FaceColor = 'r';
        flag = true;
    end

    if T(i) > t_trans(2)
        p.String = 'Back Thrust';
    else
        p.String = 'Free Fall';
    end
    drawnow
end

%% Calculation
syms x0 v0 Tmax tf m g V
t = tf - V/g;
t_trans = m*v0/Tmax + tf*(1-m*g/Tmax);
eqn = x0 + v0 * t_trans - 1/2 * g * t_trans^2 + (g - Tmax/m) * tf * (t - t_trans) - 1/2 * (g - Tmax/m)*(t^2 - t_trans^2);

eqn2 = simplify(eqn);
coeff = collect(eqn2,V);

%% Function for state propagation
function f = func(X,t,t_trans,m,Tmax)
    g = 9.81;
    if t < t_trans
        f = [X(2);-g];
    else
        f = [X(2);Tmax/m - g];
    end
end

