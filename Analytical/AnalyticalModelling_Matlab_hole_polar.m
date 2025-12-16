clear; clc; close all;

%% Parameters are defined

F = 4000;    % [N] Applied load
t = 8;       % [mm] Thickness of the specimen
d = 12;      % [mm] Diameter of hole
W = 25;      % [mm] Width of the specimen

a = d/2;               % [mm] Hole radius 
Sigma_Inf = F/(t*W);   % [MPa] Stress at infinity

%% Stress Concentration Factor is calculated
Kt = 2 + 0.284*(1 - (d/W)) - 0.600*(1 - (d/W))^2 + 1.32*(1 - (d/W))^3;  

%% Maximum stress is calculated
Sigma_max = (Kt*1/(1-d/W))*Sigma_Inf;

%% Radial stress is calculated
Sigma_rr = @(r,theta) (Sigma_Inf/2).*(1-(a./r).^2)+(Sigma_Inf/2).*(1-4*(a./r).^2+3.*(a./r).^4).*cos(2.*theta);

%% Circumferential stress is calculated
Sigma_oo = @(r,theta) (Sigma_Inf/2).*(1+(a./r).^2)-(Sigma_Inf/2).*(1+3*(a./r).^4).*cos(2.*theta);

%% Shear stress is calculated
tau_ro = @(r,theta) -(Sigma_Inf/2).*(1+2.*(a./r).^2-3.*(a./r).^4).*sin(2.*theta);

%% Von Mises stress is calculated
Sigma_von = @(r,theta) sqrt(Sigma_rr(r,theta).^2+Sigma_oo(r,theta).^2-(Sigma_rr(r,theta).*Sigma_oo(r,theta))+3*tau_ro(r,theta).^2);

%% Range for r and theta are defined
rmax = sqrt(2)*(W/2);              % reach plate corners
r = linspace(a, rmax, 400);        % r from a to 12 mm
theta = linspace(-pi, pi, 800);    % theta from 0 to 2*pi
[R, Theta] = meshgrid(r, theta);   % Create a meshgrid for r and theta
[X,Y] = pol2cart(Theta,R);         % Convert to cartesian coordinats
Theta_load = Theta - pi/2;         % Rotates the polar coordinate system so theta starts on the y axis. 
                                   % (Matlab thinks theta starts at x not y)

%% Evaluate 
Sigmar = Sigma_rr(R,Theta_load);
Sigmao = Sigma_oo(R,Theta_load);
tau_rth = tau_ro(R,Theta_load);
Sigmav = Sigma_von(R,Theta_load);

%% Stress transformation to cartesian coordinates
sigma_x = Sigmar .* cos(Theta).^2 + Sigmao .* sin(Theta).^2 - tau_rth .* sin(2*Theta);
sigma_y = Sigmar .* sin(Theta).^2 + Sigmao .* cos(Theta).^2 + tau_rth .* sin(2*Theta);
tau_xy  = (Sigmar - Sigmao) .* sin(Theta) .* cos(Theta) + tau_rth .* cos(2*Theta);

%% Calculate maximum values of stresses
max_Sigmar = max(Sigmar(:));
max_Sigmao = max(Sigmao(:));
max_tau = max(abs(tau_rth(:)));
max_Sigmav = max(Sigmav(:));

max_sigma_x = max(sigma_x(:));
max_sigma_y = max(sigma_y(:));
max_tau_xy = max(tau_xy(:));

%% Display maximum values
fprintf('Maximum Radial Stress: %.2f MPa\n', max_Sigmar);
fprintf('Maximum Circumferential Stress: %.2f MPa\n', max_Sigmao);
fprintf('Maximum Shear Stress: %.2f MPa\n', max_tau);
fprintf('Maximum von Mises Stress: %.2f MPa\n', max_Sigmav);

fprintf('Maximum stress in X direction: %.2f MPa\n', max_sigma_x);
fprintf('Maximum stress in Y direction: %.2f MPa\n', max_sigma_y);
fprintf('Maximum shear stress in X Y: %.2f MPa\n', max_tau_xy);

%% Plot Radial stresses (MPa)
figure('Color','w')
contourf(X, Y, Sigmar, 100, "LineStyle","none")
set(gcf, 'Name', 'sigma_r')
axis equal tight
xlim([-W/2 W/2]); ylim([-W/2 W/2])
xlabel('$x$ [mm]', 'FontSize', 18,'Interpreter','latex')
ylabel('$y$ [mm]', 'FontSize', 18,'Interpreter','latex')
ax = gca;
ax.FontSize = 14;
colormap("jet")
c = colorbar("AxisLocation","out");
c.Label.String = "[MPa]";
c.Label.FontSize = 14;

%% Plot Circumferential stresses (MPa)
figure('Color','w')
contourf(X, Y, Sigmao, 100, "LineStyle","none")
set(gcf, 'Name', 'sigma_th')
axis equal tight
xlim([-W/2 W/2]); ylim([-W/2 W/2])
xlabel('$x$ [mm]', 'FontSize', 18,'Interpreter','latex')
ylabel('$y$ [mm]', 'FontSize', 18,'Interpreter','latex')
ax = gca;
ax.FontSize = 14;
colormap("jet")
c = colorbar("AxisLocation","out");
c.Label.String = "[MPa]";
c.Label.FontSize = 14;

%% Plot Shear stresses (MPa)
figure('Color','w')
contourf(X, Y, tau_rth, 100, "LineStyle","none")
set(gcf, 'Name', 'tau_rth')
axis equal tight
xlim([-W/2 W/2]); ylim([-W/2 W/2])
xlabel('$x$ [mm]', 'FontSize', 18,'Interpreter','latex')
ylabel('$y$ [mm]', 'FontSize', 18,'Interpreter','latex')
ax = gca;
ax.FontSize = 14;
colormap("jet")
c = colorbar("AxisLocation","out");
c.Label.String = '[MPa]';
c.Label.FontSize = 14;

%% Plot von Mises stresses (MPa)
figure('Color','w')
contourf(X, Y, Sigmav, 100, 'LineStyle', 'none')
set(gcf, 'Name', 'sigma_von')
axis equal tight
xlim([-W/2 W/2]); ylim([-W/2 W/2])
xlabel('$x$ [mm]', 'FontSize', 18,'Interpreter','latex')
ylabel('$y$ [mm]', 'FontSize', 18,'Interpreter','latex')
ax = gca;
ax.FontSize = 14;
colormap("jet")
c = colorbar("AxisLocation","out");
c.Label.String = '[MPa]';
c.Label.FontSize = 14;

%% Plot x stresses (MPa)
figure('Color','w')
contourf(X, Y, sigma_x, 100, "LineStyle", "none")
set(gcf, 'Name', 'sigma_x')
axis equal tight
xlim([-W/2 W/2]); ylim([-W/2 W/2])
xlabel('$x$ [mm]', 'FontSize', 18,'Interpreter','latex')
ylabel('$y$ [mm]', 'FontSize', 18,'Interpreter','latex')
ax = gca;
ax.FontSize = 14;
colormap("jet")
c = colorbar("AxisLocation","out");
c.Label.String = '[MPa]';
c.Label.FontSize = 14;

%% Plot y stresses (MPa)
figure('Color','w')
contourf(X, Y, sigma_y, 100, "LineStyle", "none")
set(gcf, 'Name', 'sigma_y')
axis equal tight
xlim([-W/2 W/2]); ylim([-W/2 W/2])
xlabel('$x$ [mm]', 'FontSize', 18,'Interpreter','latex')
ylabel('$y$ [mm]', 'FontSize', 18,'Interpreter','latex')
ax = gca;
ax.FontSize = 14;
colormap("jet")
c = colorbar("AxisLocation","out");
c.Label.String = '[MPa]';
c.Label.FontSize = 14;

%% Plot shear stress (MPa)
figure('Color','w')
contourf(X, Y, tau_xy, 100, "LineStyle", "none")
set(gcf, 'Name', 'tau_xy')
axis equal tight
xlim([-W/2 W/2]); ylim([-W/2 W/2])
xlabel('$x$ [mm]', 'FontSize', 18,'Interpreter','latex')
ylabel('$y$ [mm]', 'FontSize', 18,'Interpreter','latex')
ax = gca;
ax.FontSize = 14;
colormap("jet")
c = colorbar("AxisLocation","out");
c.Label.String = '[MPa]';
c.Label.FontSize = 14;

%% Line for comparison plots
% Start and end points are defined in cartesian coordinates
x_start = -12.5;    x_end = 12.5;
y_start = 8.0;      y_end = 8.0;

% Difference between start and end points 
dx = x_end-x_start;
dy = y_end-y_start;

% Line parameter to make N points along the line
N = 200;
t = linspace(0,1,N);

% Calculate the stress values along the defined line
x_line = x_start + t * dx;
y_line = y_start + t * dy;

% Calculate distance along line
s_line = sqrt(dx^2+dy^2)*t;

% Interpolation between all points.
F_sig_x = scatteredInterpolant(X(:), Y(:), sigma_x(:), 'linear', 'none');
F_sig_y = scatteredInterpolant(X(:), Y(:), sigma_y(:), 'linear', 'none');
F_tau_xy = scatteredInterpolant(X(:), Y(:), tau_xy(:), 'linear', 'none');
F_sig_eq = scatteredInterpolant(X(:), Y(:), Sigmav(:), 'linear', 'none');

% Stresses on the line are calculated
sigma_x_line = F_sig_x(x_line,y_line);
sigma_y_line = F_sig_y(x_line,y_line);
tau_xy_line = F_tau_xy(x_line,y_line);
sigma_eq_line = F_sig_eq(x_line,y_line);


% Compile data and export
Stress_lines = { ...
    sigma_x_line, ...
    sigma_y_line, ...
    tau_xy_line, ...
    sigma_eq_line, ...
    };

names = {'sigma_x','sigma_y','tau_xy','sigma_eq'};
outDir = 'Lineout_data';


for k = 1:numel(Stress_lines)
    
    T{k} = table(s_line(:), Stress_lines{k}(:), 'VariableNames', {'Length_mm', names{k}});

    filename = sprintf('Hole_Kirsch_%s.csv', names{k});
    filepath = fullfile(outDir, filename);

    
    writetable(T{k}, filepath);

    fprintf('Saved lineout data as: %s\n', filepath);
end
