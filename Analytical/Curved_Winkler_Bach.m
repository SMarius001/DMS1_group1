clear; clc; close all
syms th_deg r real

% Parameters are defined
P = 4000;       % [N] Applied load
Rct = 45e-3;      % [m] Radius to center
b = 8e-3;    % [m] Thickness of the specimen
h = 25e-3;      % [m] Distance from inner to outer edge 

A = (8e-3)*h; % [m^2] Cross sectional area
I = 1/12 * h^3*b;  % Moment of inertia

%% Inner forces are calculated
N(th_deg) = P*sind(th_deg);     % Normal force 
V(th_deg) = P*cosd(th_deg);     % Shear force
M(th_deg) = Rct*sind(th_deg)*(P);        % Moment


%% The Shear stress is calculated
y(r) = r-Rct;
tau_xy(th_deg,r) = V(th_deg)/(2*I)*(h^2/4-y(r)^2);


%% The tangential stress is calculated using Winkler Bachs formula
ri = Rct - h/2;   % [m] Inner radius 
ro = Rct + h/2;   % [m] Outer radius
R = h/log(ro/ri);   % neutral axis radius
e = R - Rct; 
% e is flipped compared to Advanced Mechanics of Materials, 
% since theta is assumed to run from 0deg at the top of the beam to 180deg at the bottom running clockwise,
% This means that the sign is flipped compared to the book where theta runs counterclockwise 

% Axial and bending stresses are calculated seperately
sig_axial(th_deg) = N(th_deg) / A;
sig_bend(th_deg, r) = -M(th_deg)*(R-r)/(A*e*r);

% Superposition is used to get the tangential stresses 
sig_th(th_deg, r) = sig_axial(th_deg) + sig_bend(th_deg, r);


%% Pricipal stress is calculated
sig_r = sym(0);

sig_p1(th_deg,r) = (sig_th+sig_r)/2 + sqrt(((sig_th-sig_r)/2)^2+tau_xy^2);
sig_p2(th_deg,r) = (sig_th+sig_r)/2 - sqrt(((sig_th-sig_r)/2)^2+tau_xy^2);


%% Von Mises stress is calculated
sig_vm(th_deg, r) = sqrt(sig_th(th_deg, r).^2 + 3*tau_xy(th_deg, r).^2);
sig_vmp(th_deg,r) = sqrt(sig_p1(th_deg,r)^2 + sig_p2(th_deg,r)^2 - sig_p1(th_deg,r)*sig_p2(th_deg,r));

%% Numeric handle
sigth_fun = matlabFunction(sig_th, 'Vars', [th_deg, r]);
tau_fun = matlabFunction(tau_xy, 'Vars', [th_deg, r]);
sigp1_fun = matlabFunction(sig_p1, 'Vars', [th_deg, r]);
sigp2_fun = matlabFunction(sig_p2, 'Vars', [th_deg, r]);
sig_vm_fun = matlabFunction(sig_vm, 'Vars', [th_deg, r]);
sig_vmp_fun = matlabFunction(sig_vmp, 'Vars', [th_deg, r]);

%%
% Grids
th_span = linspace(0,180,721);
r_span  = linspace(ri, ro, 241);
[THdeg, Rgrid] = meshgrid(th_span, r_span);

% Evaluate (Pa)
Sigma_th = sigth_fun(THdeg, Rgrid);
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

% Plot rangential stresses (MPa)
figure('Color','w');
contourf(Xmm, Ymm, Sigma_th/1e6, 100, 'LineColor','none');
axis equal tight
xlabel('{\it x} [mm]'); ylabel('{\it y} [mm]');
cb = colorbar; ylabel(cb,'[MPa]'); 
% title('Tangential stresses \sigma_\theta – Winkler–Bach')
set(gcf, 'Name', 'WB_Tangential');   % used for filename if no title
set(gca, 'Tag',  'tg_WB');           % secondary fallback for filename
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
%title('Principal stress \sigma_1 – Winkler Bach')
set(gcf, 'Name', 'WB_Principal1');   % used for filename if no title
set(gca, 'Tag',  'p1_WB');           % secondary fallback for filename
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
%title('Principal stress \sigma_2 – Winkler Bach')
set(gcf, 'Name', 'WB_Principal2');   % used for filename if no title
set(gca, 'Tag',  'p2_WB');           % secondary fallback for filename
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
% title('Shear stresses \tau_{xy} – Wikler Bach')
set(gcf, 'Name', 'WB_Shear');   % used for filename if no title
set(gca, 'Tag',  'sh_WB');           % secondary fallback for filename
hold on
thp = linspace(min(th_span), max(th_span), 400);
plot(1e3*ri*sind(thp), 1e3*ri*cosd(thp), 'k-', 'LineWidth', 1.2);
plot(1e3*ro*sind(thp), 1e3*ro*cosd(thp), 'k-', 'LineWidth', 1.2);
colormap(jet)
hold off

% Plot von Mieses stresses (MPa)
figure('Color','w');
contourf(Xmm, Ymm, Sigma_vm/1e6, 100, 'LineColor','none');
axis equal tight
xlabel('{\it x} [mm]'); ylabel('{\it y} [mm]');
cb = colorbar; ylabel(cb,'[MPa]');
% title('von Mises stresses \sigma_{vm} – Winkler Bach')
set(gcf, 'Name', 'WB_von_Mises');   % used for filename if no title
set(gca, 'Tag',  'vm_WB');           % secondary fallback for filename
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
%title('Principal von Mises stress \sigma_{vm} – Winkler Bach')
set(gcf, 'Name', 'WB_von_Mises_p');   % used for filename if no title
set(gca, 'Tag',  'vm_WBp');          % secondary fallback for filename
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
RES.method   = 'Winkler-Bach';          
RES.units    = 'Pa';
RES.THdeg    = THdeg;                   % theta grid (deg), size [nr x nt]
RES.Rgrid    = Rgrid;                   % radius grid (m), size [nr x nt]
RES.Xmm      = Xmm;                     % plotting grid (mm) (optional)
RES.Ymm      = Ymm;

RES.sigma_th = Sigma_th;                % hoop/tangential normal stress [Pa]
RES.tau      = Tau_xy;                  % in-plane shear τ_{rθ} [Pa]
RES.sigma_p1 = Sigma_p1;                % Principal stress 1 [Pa]
RES.sigma_p2 = Sigma_p2;                % Principal stress 2 [Pa]
RES.sigma_vm = Sigma_vm;                % von Mises (plane stress) [Pa]

RES.ri       = ri;                      % [m]
RES.ro       = ro;                      % [m]
RES.Rc       = Rct;                       % centroid radius [m] (if you have it)

save('WB_result.mat','RES');


