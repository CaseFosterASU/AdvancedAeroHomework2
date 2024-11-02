% Clear previous variables and close figures
clear; close all; clc;

% Define constants
U_free = 1;   % Freestream velocity
chord = 1;    % Chord length
eps = 0.01;   % Small offset for flatness
r = chord/2;  % Radius of circular cylinder

% Create circular points (f-plane)
num_points = 200;               % Number of points around the circle
theta_vals = linspace(0, 2*pi, num_points);
circle_points = r * exp(1i * theta_vals);

% Apply Joukowski transformation (z-plane)
transformed_points = circle_points + (r^2)./circle_points;

% Center the plate along the x-axis
x_shift = (max(real(transformed_points)) + min(real(transformed_points))) / 2;
transformed_points = transformed_points - x_shift;

% Stagnation points
stagnation_points = [r, -r];

% Plot the circle in the f-plane
figure;
plot(real(circle_points), imag(circle_points), 'b-', 'LineWidth', 1.5); % Blue circle
hold on;
plot(real(stagnation_points), imag(stagnation_points), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8); % Green points
plot([-2, 2], [0, 0], 'k--'); % x-axis
plot([0, 0], [-2, 2], 'k--'); % y-axis
axis equal;
xlabel('Real Axis');
ylabel('Imaginary Axis');
title('Circle (f-plane)');
xlim([-2, 2]);
ylim([-2, 2]);

% Plot the transformed shape in the z-plane
figure;
plot(real(transformed_points), imag(transformed_points), 'r-', 'LineWidth', 1.5); % Red line for flat plate
hold on;
plot([-2, 2], [0, 0], 'k--'); % x-axis
plot([0, 0], [-2, 2], 'k--'); % y-axis
axis equal;
xlabel('Real Axis');
ylabel('Imaginary Axis');
title('Flat Plate (z-plane)');
xlim([-2, 2]);
ylim([-2, 2]);

% Basic annotation for verification
text(-1.5, 1.5, ['Center of plate: x = ', num2str(mean(real(transformed_points)))], 'FontSize', 8);


%% Circular arc
u_inf = 1;                  % Free stream velocity (not used in this transformation)
beta = linspace(0, pi, 100); % Angle for parameterizing the circle
c = 1;                      % Chord length of the airfoil
a = c / 2;                  % Radius of the circle (based on chord length)
m = 1;                % Camber offset (adjust this for more or less camber)
stagnation_points = [a, -a];

% Define the circle in the f-plane with a camber offset
mu = i * m;                % Vertical offset for camber
f = a * exp(1i * beta) + mu; % Circle in the f-plane, offset vertically by mu
f_low = a * exp(1i * beta) - mu;
% Apply the Joukowski transformation to create the circular arc airfoil
f_z = f + (c^2) ./ f;

% Plot the original circle in the f-plane
figure;
hold on;
plot(real(f), imag(f), 'b-', 'LineWidth', 1.5);
plot(real(f_low), -imag(f_low), 'b-', 'LineWidth', 1.5);
plot(real(mu), imag(mu), 'ro', 'MarkerFaceColor', 'r'); % Mark the center of the circle
xlabel('Real(f)');
plot(real(stagnation_points), imag(stagnation_points), 'go', ...
    'MarkerFaceColor', 'g', 'MarkerSize', 10);
ylabel('Imaginary(f)');
title('Circle in the f-plane (for Circular Arc)');
axis equal;
grid on;

% Plot the transformed airfoil in the z-plane
figure;
plot(real(f_z), imag(f_z), 'r-', 'LineWidth', 1.5);
xlabel('Real(z)');
ylabel('Imaginary(z)');
title('Circular Arc Airfoil in the z-plane');
axis equal;
grid on;
%% Joukouski Symmetric Airfoil 
% Joukowski Symmetric Airfoil
c = 1;
epsi = 0.1;
theta = linspace(0, 2*pi, 1000);
a = (c / 4) * (1 + epsi);
mu = -epsi * c / 4;

% Define Circle and Apply Joukowski Transform
f = a * exp(1i * theta) + mu;
f_z = f + c^2 ./ (16 * f);

stagnation_points = [a, -a];

% Plot Circle Plane
figure;
plot(real(f), imag(f), 'b-', 'LineWidth', 1.5);
hold on;
plot(real(stagnation_points), imag(stagnation_points), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
title('Circle Plane');
xlabel('Real Axis');
ylabel('Imaginary Axis');
axis equal;
grid on;

% Plot Transformed Airfoil Plane
figure;
plot(real(f_z), imag(f_z), 'r-', 'LineWidth', 1.5);
title('Transformed Airfoil Plane');
xlabel('Real Axis');
ylabel('Imaginary Axis');
axis equal;
grid on;

