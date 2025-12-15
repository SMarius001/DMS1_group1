clear; clc; close all


syms r th_deg
Pf = -4000;  % [N] Applied point load (Negative since it points the opposite way in Theory of Elasticity)
h = 25e-3;   % [m] Distance from inner to outer edge
R = 45e-3;   % [m] Radius to the centroid
t = 8e-3;    % [m] Thickness of the specimen

ri = R - h/2;   % [m] Inner radius 
ro = R + h/2;   % [m] Outer radius
a = ri; % Converting to the Theory of Elasticity notation.
b = ro; % Converting to the Theory of Elasticity notation.
P = Pf/t; % [N/m] The point load is converted to a load pr unit thickness 


%% Constants of the solution are defined
N = a^2 - b^2 + (a^2 + b^2)*log(b/a);

A = P/(2*N);
B = - P*a^2*b^2/(2*N);
D = - P/N*(a^2+b^2);

%% Stress components are calculated
sig_r = (2*A*r-2*B/r^3+D/r)*sind(th_deg);
sig_th = (6*A*r+2*B/r^3+D/r)*sind(th_deg);
tau_rth = -(2*A*r-2*B/r^3+D/r)*cosd(th_deg);

%% Pricipal stress is calculated
sig_p1(th_deg,r) = (sig_th+sig_r)/2 + sqrt(((sig_th-sig_r)/2)^2+tau_rth^2);
sig_p2(th_deg,r) = (sig_th+sig_r)/2 - sqrt(((sig_th-sig_r)/2)^2+tau_rth^2);

%% von Mises stress is calculated
sig_vm = sqrt(sig_th^2+sig_r^2-sig_th*sig_r+3*tau_rth^2);
sig_vmp(th_deg,r) = sqrt(sig_p1(th_deg,r)^2 + sig_p2(th_deg,r)^2 - sig_p1(th_deg,r)*sig_p2(th_deg,r));


%% Numeric handle
sig_fun = matlabFunction(sig_th, 'Vars', [th_deg, r]);
sig_r_fun = matlabFunction(sig_r, 'Vars', [th_deg,r]);
tau_fun = matlabFunction(tau_rth, 'Vars', [th_deg, r]);
sigp1_fun = matlabFunction(sig_p1, 'Vars', [th_deg, r]);
sigp2_fun = matlabFunction(sig_p2, 'Vars', [th_deg, r]);
sig_vm_fun = matlabFunction(sig_vm, 'Vars', [th_deg, r]);
sig_vmp_fun = matlabFunction(sig_vmp, 'Vars', [th_deg, r]);

%% Plotting
% Grids
th_span = linspace(0,180,180);
r_span  = linspace(ri, ro, 100);
[THdeg, Rgrid] = meshgrid(th_span, r_span);

% Evaluate (Pa)
Sigma_th = sig_fun(THdeg, Rgrid);
Sigma_r = sig_r_fun(THdeg, Rgrid);
Tau_xy = tau_fun(THdeg,Rgrid); 
Sigma_p1 = sigp1_fun(THdeg, Rgrid);
Sigma_p2 = sigp2_fun(THdeg, Rgrid);
Sigma_vm = sig_vm_fun(THdeg,Rgrid);
Sigma_vmp = sig_vmp_fun(THdeg,Rgrid);


% Rotate 90° clockwise for your preferred view
Xmm = 1e3 * (Rgrid .* sind(THdeg));
Ymm = 1e3 * (Rgrid .* cosd(THdeg));

% ---- mark max & min on the map ----
[SigMax, iMax] = max(Sigma_vm(:));
[SigMaxP, iMaxP] = max(Sigma_vmp(:));

xMax = Xmm(iMax);  yMax = Ymm(iMax);
xMaxP = Xmm(iMaxP);  yMaxP = Ymm(iMaxP);

% Plot Tangential stresses (MPa)
figure('Color','w');
contourf(Xmm, Ymm, Sigma_th/1e6, 100, 'LineColor','none');
axis equal tight
xlabel('{\it x} [mm]'); ylabel('{\it y} [mm]');
cb = colorbar; ylabel(cb,'[MPa]'); 
%title('Tangential stresses \sigma_\theta – Airys')
set(gcf, 'Name', 'Air_Tangential');   % used for filename if no title
set(gca, 'Tag',  'tg_Air');           % secondary fallback for filename
hold on
thp = linspace(min(th_span), max(th_span), 400);
plot(1e3*ri*sind(thp), 1e3*ri*cosd(thp), 'k-', 'LineWidth', 1.2);
plot(1e3*ro*sind(thp), 1e3*ro*cosd(thp), 'k-', 'LineWidth', 1.2);
colormap(jet)
hold off

% Plot Radial stresses (MPa)
figure('Color','w');
contourf(Xmm, Ymm, Sigma_r/1e6, 100, 'LineColor','none');
axis equal tight
xlabel('{\it x} [mm]'); ylabel('{\it y} [mm]');
cb = colorbar; ylabel(cb,'[MPa]'); 
%title('Radial stresses \sigma_r – Airys')
set(gcf, 'Name', 'Air_Radial');   % used for filename if no title
set(gca, 'Tag',  'rd_Air');           % secondary fallback for filename
hold on
thp = linspace(min(th_span), max(th_span), 400);
plot(1e3*ri*sind(thp), 1e3*ri*cosd(thp), 'k-', 'LineWidth', 1.2);
plot(1e3*ro*sind(thp), 1e3*ro*cosd(thp), 'k-', 'LineWidth', 1.2);
colormap(jet)
hold off

% Plot Pricipal stress 1
figure('Color','w');
contourf(Xmm, Ymm, Sigma_p1/1e6, 100, 'LineColor','none');
axis equal tight
xlabel('{\it x} [mm]'); ylabel('{\it y} [mm]');
cb = colorbar; ylabel(cb,'[MPa]'); 
%title('Principal stress \sigma_1 – Airys')
set(gcf, 'Name', 'Air_Principal1');   % used for filename if no title
set(gca, 'Tag',  'p1_Air');           % secondary fallback for filename
hold on
thp = linspace(min(th_span), max(th_span), 400);
plot(1e3*ri*sind(thp), 1e3*ri*cosd(thp), 'k-', 'LineWidth', 1.2);
plot(1e3*ro*sind(thp), 1e3*ro*cosd(thp), 'k-', 'LineWidth', 1.2);
colormap(jet)
hold off

% Plot Pricipal stress 2
figure('Color','w');
contourf(Xmm, Ymm, Sigma_p2/1e6, 100, 'LineColor','none');
axis equal tight
xlabel('{\it x} [mm]'); ylabel('{\it y} [mm]');
cb = colorbar; ylabel(cb,'[MPa]'); 
%title('Principal stress \sigma_2 – Airys')
set(gcf, 'Name', 'Air_Principal2');   % used for filename if no title
set(gca, 'Tag',  'p2_Air');           % secondary fallback for filename
hold on
thp = linspace(min(th_span), max(th_span), 400);
plot(1e3*ri*sind(thp), 1e3*ri*cosd(thp), 'k-', 'LineWidth', 1.2);
plot(1e3*ro*sind(thp), 1e3*ro*cosd(thp), 'k-', 'LineWidth', 1.2);
colormap(jet)
hold off

% Plot Shear stresses (MPa)
figure('Color','w');
contourf(Xmm, Ymm, Tau_xy/1e6, 100, 'LineColor','none');
axis equal tight
xlabel('{\it x} [mm]'); ylabel('{\it y} [mm]');
cb = colorbar; ylabel(cb,'[MPa]'); 
%title('Shear stresses \tau_{xy} – Airys')
set(gcf, 'Name', 'Air_Shear');   % used for filename if no title
set(gca, 'Tag',  'sh_Air');           % secondary fallback for filename
hold on
thp = linspace(min(th_span), max(th_span), 400);
plot(1e3*ri*sind(thp), 1e3*ri*cosd(thp), 'k-', 'LineWidth', 1.2);
plot(1e3*ro*sind(thp), 1e3*ro*cosd(thp), 'k-', 'LineWidth', 1.2);
colormap(jet)
hold off

figure('Color','w');
contourf(Xmm, Ymm, Sigma_vm/1e6, 100, 'LineColor','none');
axis equal tight
xlabel('{\it x} [mm]'); ylabel('{\it y} [mm]');
cb = colorbar; ylabel(cb,'[MPa]');
%title('von Mises stresses \sigma_{vm} – Airys')
set(gcf, 'Name', 'Air_von_Mises');   % used for filename if no title
set(gca, 'Tag',  'vm_Air');           % secondary fallback for filename
hold on
thp = linspace(min(th_span), max(th_span), 400);
colormap(jet)
plot(1e3*ri*sind(thp), 1e3*ri*cosd(thp), 'k-', 'LineWidth', 1.2);
plot(1e3*ro*sind(thp), 1e3*ro*cosd(thp), 'k-', 'LineWidth', 1.2);

% peak marker + label
plot(xMax, yMax, 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor','w','MarkerEdgeColor','k');
dx = 2; % mm left offset
hTxt = text(xMax - dx, yMax, sprintf('%.1f MPa', SigMax/1e6), ...
    'HorizontalAlignment','right','VerticalAlignment','middle', ...
    'BackgroundColor','w','Margin',2,'FontWeight','bold', ...
    'Tag','peakLabel', ...   % <-- tag so the saver can find it
    'FontSize',14);          % base size if you run without the saver
hold off

% Plot principal von Mises stress
figure('Color','w');
contourf(Xmm, Ymm, Sigma_vmp/1e6, 100, 'LineColor','none');
axis equal tight
xlabel('{\it x} [mm]'); ylabel('{\it y} [mm]');
cb = colorbar; ylabel(cb,'[MPa]');
%title('Principal von Mises stress \sigma_{vm} – Airys')
set(gcf, 'Name', 'Air_von_Mises_p');   % used for filename if no title
set(gca, 'Tag',  'vm_Airp');          % secondary fallback for filename
hold on
thp = linspace(min(th_span), max(th_span), 400);
colormap(jet)
plot(1e3*ri*sind(thp), 1e3*ri*cosd(thp), 'k-', 'LineWidth', 1.2);
plot(1e3*ro*sind(thp), 1e3*ro*cosd(thp), 'k-', 'LineWidth', 1.2);
% peak marker + label
plot(xMaxP, yMaxP, 'o', 'MarkerSize', 5, ...
    'MarkerFaceColor','w','MarkerEdgeColor','k');
dx = 2; % mm left offset
hTxtp = text(xMaxP - dx, yMaxP, sprintf('%.1f MPa', SigMaxP/1e6), ...
    'HorizontalAlignment','right','VerticalAlignment','middle', ...
    'BackgroundColor','w','Margin',2,'FontWeight','bold', ...
    'Tag','peakLabel', ...   % <-- tag so the saver can find it
    'FontSize',14);          % base size if you run without the saver
hold off

%% Export results for comparison script
RES.method   = 'Airys';          
RES.units    = 'Pa';
RES.THdeg    = THdeg;                   % theta grid (deg), size [nr x nt]
RES.Rgrid    = Rgrid;                   % radius grid (m), size [nr x nt]
RES.Xmm      = Xmm;                     % plotting grid (mm) (optional)
RES.Ymm      = Ymm;

RES.sigma_th = Sigma_th;                % hoop/tangential normal stress [Pa]
RES.tau      = Tau_xy;                  % in-plane shear τ_{rθ} [Pa]
RES.sigma_r  = Sigma_r;                 % Radial normal stress [Pa]
RES.sigma_p1 = Sigma_p1;                % Principal stress 1 [Pa]
RES.sigma_p2 = Sigma_p2;                % Principal stress 2 [Pa]
RES.sigma_vm = Sigma_vm;                % von Mises (plane stress) [Pa]

RES.ri       = ri;                      % [m]
RES.ro       = ro;                      % [m]
RES.Rc       = R;                       % centroid radius [m] (if you have it)

save('Air_result.mat','RES');

