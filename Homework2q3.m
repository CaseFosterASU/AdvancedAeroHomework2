%% Question 3 
clc; clear;

% Geometry and Flow Parameters
c = 1;
epsi = 0.1;
theta = linspace(0, 2*pi, 128); % Angular positions in radians
a = c/4 * (1 + epsi);           % Radius for symmetric airfoil
beta = 0;                       % Camber angle (radians)
aoa = 0;                        % Angle of attack (radians)
mue = -epsi * c / 4;            % Offset for symmetric airfoil
u_inf = 1;                      % Freestream velocity
gamma = 4 * pi * a * u_inf * sin(aoa + beta); % Circulation (Gamma)

% Circle and Airfoil Geometry in Complex Planes
f = a * exp(1i * theta) + mue;  % Circle in f-plane
f_z = f + c^2 ./ (16 * f);      % Transformed airfoil in z-plane
x = (c / 2) * cos(theta);       % x-coordinates along the airfoil

% Streamlines in Circle Plane (f-plane)
x_r = [-2, 2];
z_r = [-2, 2];
[X, Z] = meshgrid(linspace(x_r(1), x_r(2), 100), linspace(z_r(1), z_r(2), 100));
x_stream = -2 * ones(size(X(:, 1))); % Starting x-position for streamlines
z_stream = Z(:, 1);                  % Starting y-positions for streamlines

% Compute Complex Potential in Circle Plane
w_circle = zeros(size(X));
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        kk = X(i, j) + 1i * Z(i, j);
        w_circle(i, j) = (u_inf * exp(-1i * aoa)) + ((1i * gamma) / (2 * pi)) * (1 / (kk - mue)) ...
                         - (u_inf * (a^2)) * (exp(1i * aoa) / (kk - mue)^2);
    end
end

% Plot Streamlines in Circle Plane (f-plane)
figure;
plot(real(f), imag(f), 'r-', 'LineWidth', 1.5); % Circle in f-plane
hold on;
streamline(X, Z, real(w_circle), -imag(w_circle), x_stream, z_stream);
title('Streamlines around Circle in f-plane at 0 Degrees Angle of Attack');
xlabel('X');
ylabel('Z');
axis equal;
grid on;
hold off;

% Compute Complex Potential for Airfoil Plane (z-plane) with Joukowski Transformation
w_airfoil = zeros(size(X));
for i = 1:size(X, 1)
    for j = 1:size(X, 2)
        ww = X(i, j) + 1i * Z(i, j);
        
        % Inverse Joukowski Transformation
        if real(ww) < 0
            kk = 0.5 * (ww - sqrt(ww^2 - (c^2 / 4)));
        else
            kk = 0.5 * (ww + sqrt(ww^2 - (c^2 / 4)));
        end
        
        % Complex potential in f-plane and transformed in z-plane
        w_f = (u_inf * exp(-1i * aoa)) + ((1i * gamma) / (2 * pi)) * (1 / (kk - mue)) ...
              - (u_inf * (a^2)) * (exp(1i * aoa) / (kk - mue)^2);
        w_airfoil(i, j) = w_f / (1 - (c^2 / (16 * kk^2))); % Adjusted potential in z-plane
    end
end

% Plot Streamlines in Airfoil Plane (z-plane)
figure;
plot(real(f_z), imag(f_z), 'r-', 'LineWidth', 1.5); % Airfoil in z-plane
hold on;
streamline(X, Z, real(w_airfoil), -imag(w_airfoil), x_stream, z_stream);
title('Streamlines around Airfoil in z-plane at 0 Degrees Angle of Attack');
xlabel('X');
ylabel('Z');
axis equal;
grid on;
hold off;

% Velocity and Pressure Coefficients on Airfoil Surface
k = exp(-1i * theta); % Complex exponential factor
W_f = 2i * u_inf .* k .* (sin(aoa + beta) - sin(aoa - theta)); % Complex potential in f-plane
W_y = W_f ./ (1 - (c^2 ./ (16 * f.^2))); % Adjusted potential for z-plane

% Velocity Magnitude and Pressure Coefficient
U = real(W_y);
W = imag(W_y);
v_t = sqrt(U.^2 + W.^2); % Velocity magnitude on airfoil surface
cp = 1 - (v_t / u_inf).^2; % Pressure coefficient

% Plot Velocity Distribution as a Function of x
figure;
plot(x, v_t, 'g-', 'LineWidth', 1.5); % Blue line for velocity distribution
xlabel('X Coordinate');
ylabel('Velocity Distribution');
title('Velocity Distribution along the Airfoil Surface');
grid on;

% Plot Pressure Coefficient Distribution as a Function of x
figure;
plot(x, -cp, 'g-', 'LineWidth', 1.5); % Green line for pressure distribution
xlabel('X Coordinate');
ylabel('Pressure Coefficient (Cp)');
title('Pressure Coefficient Distribution along the Airfoil Surface');
grid on;
