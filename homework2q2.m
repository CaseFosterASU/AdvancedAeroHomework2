% Clear workspace and close figures
clear all; close all; clc;

% Define parameters
c = 1;                   % Chord length
beta = 2 * pi / 180;     % Camber angle
alpha = 4 * pi / 180;    % Angle of attack
u_inf = 1;               % Freestream velocity
a = (c / 4) * sec(beta); % Radius
m = a * sin(beta);       % Camber offset
mu = 1i * m;             % Center of circle
gamma = 4 * pi * a * u_inf * sin(alpha + beta); % Circulation

% Define grid for flow field
x_range = [-2, 2];
z_range = [-2, 2];
[X, Z] = meshgrid(linspace(x_range(1), x_range(2), 100), linspace(z_range(1), z_range(2), 100));

% Initialize complex potential arrays for circle plane and airfoil plane
w_circle = zeros(size(X));
w_airfoil = zeros(size(X));

% Compute complex potential for flow around the circle (f-plane)
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        ff = X(i, j) + 1i * Z(i, j);
        w_circle(i, j) = (u_inf * exp(-1i * alpha)) + ((1i * gamma) / (2 * pi)) * (1 / (ff - mu)) ...
                         - (u_inf * (a^2)) * (exp(1i * alpha) / (ff - mu)^2);
    end
end

% Streamline plot for the circle (f-plane)
figure();
streamline(X, Z, real(w_circle), -imag(w_circle), ones(size(X(:, 1)))-2, Z(:, 1));
title('Streamlines in the Circle Plane (f-plane) at 2 Degrees');
xlabel('X');
ylabel('Z');
axis equal;
grid on;

% Compute complex potential for the airfoil plane (z-plane) using the Joukowski transformation
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        tt = X(i, j) + 1i * Z(i, j);
        
        % Joukowski inverse transformation to find points in the f-plane
        if real(tt) < 0
            ff = 0.5 * (tt - sqrt(tt^2 - (c^2) / 4));
        else
            ff = 0.5 * (tt + sqrt(tt^2 - (c^2) / 4));
        end
        
        % Complex potential in f-plane and adjusted in z-plane
        w_airfoil_f = (u_inf * exp(-1i * alpha)) + ((1i * gamma) / (2 * pi)) * (1 / (ff - mu)) ...
                      - (u_inf * (a^2)) * (exp(1i * alpha) / (ff - mu)^2);
        w_airfoil(i, j) = w_airfoil_f / (1 - (c^2 / (16 * ff^2))); % Adjusted velocity in z-plane
    end
end

% Streamline plot for the airfoil (z-plane)
figure();
streamline(X, Z, real(w_airfoil), -imag(w_airfoil), ones(size(X(:, 1)))-2, Z(:, 1));
title('Streamlines in the Airfoil Plane (z-plane) at 2 Degrees');
xlabel('X');
ylabel('Z');
axis equal;
grid on;
