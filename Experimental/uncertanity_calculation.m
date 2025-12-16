% Strain & Stress Uncertainty
% Single Gauge, Angle (tee-rosette), Rectangular Rosette (0/45/90)
% ================================================================
clear; clc; close all;

%% 1) MATERIAL CONSTANTS
E = 205e3; % Young's modulus [MPa]
sE = 5.31; % Absolute uncertainty in E [MPa]
nu0 = 0.285; % Poisson (calibration object) for corrections
nu_spec = 0.253; % Poisson (specimen) for stress
s_nu = 0.016; % Absolute uncertainty in specimen Poisson

%% 2) GAUGE FACTORS AND TRANSVERSE SENSITIVITIES
% --- Single Gauge ---
k_sg1 = 2.06; % Gauge factor
k_sg2 = 2.06;
s_k_sg1_rel = 0.0; % tolerance
s_k_sg2_rel = 0.0;
kt_sg1 = 0.008; % 0.8%transverse sensitivity
kt_sg2 = 0.008;
s_kt_sg1_abs = 0.002; % ±0.2 percentage points
s_kt_sg2_abs = 0.002;

% --- Angle Gauge (tee-rosette) ---
k_ang1 = 1.960; % Grid 1 gauge factor
s_k_ang1_rel = 0.005; % ±0.5%
k_ang2 = 1.990; % Grid 2 gauge factor
s_k_ang2_rel = 0.005; % ±0.5%
kt_ang = 0.026; % 2.6% transverse sensitivity
s_kt_ang_abs = 0.002; % ±0.2 percentage points

% --- Rectangular Rosette (0/45/90) ---
k_ros1 = 2.09; % All three grids have same
k_ros2 = 2.09;
k_ros3 = 2.09;
s_k_ros_rel = 0.01; % ±1%
kt_ros = 0.0; % 0.0%transverse sensitivity
s_kt_ros_abs = 0.0; % since kt_ros = 0

%% 3) MEASURED STRAINS (IN STRAIN)
% Single gauge (uniaxial) strain
R_sg1 = 1419.397e-6;
R_sg2 = 55.387e-6;

% Angle gauge (two orthogonal grids)
R_ang_x = -37.805e-6; % 90° grid strain
R_ang_y = 341.242e-6; % 0° grid strain

% Rosette (0°, 45°, 90° raw strains)
R_ros_0 = 214.507e-6;
R_ros_45 = 4.988e-6;
R_ros_90 = -114.664e-6;

%% 4) RAW STRAIN UNCERTAINTIES (from k tolerances)
% u(ε_raw) = ε_raw * s_k_rel
u_sg1_raw    = abs(R_sg1)    * s_k_sg1_rel;
u_sg2_raw    = abs(R_sg2)    * s_k_sg2_rel;
u_ang_xraw   = abs(R_ang_x)  * s_k_ang1_rel; % grid 1 uses k_ang1 tolerance
u_ang_yraw   = abs(R_ang_y)  * s_k_ang2_rel; % grid 2 uses k_ang2 tolerance
u_ros_0raw   = abs(R_ros_0)  * s_k_ros_rel;  % grid 1
u_ros_45raw  = abs(R_ros_45) * s_k_ros_rel;  % grid 2
u_ros_90raw  = abs(R_ros_90) * s_k_ros_rel;  % grid 3

%% 5) SINGLE GAUGE — corrected strain and uncertainty
% ε_l = ε_raw * (1 - nu0*kt)/(1 - nu_spec*kt)
eps_sg1 = R_sg1 * (1 - nu0*kt_sg1) / (1 - nu_spec*kt_sg1);
eps_sg2 = R_sg2 * (1 - nu0*kt_sg2) / (1 - nu_spec*kt_sg2);

% Sensitivities for single gauge strain: Eq = R * (1 - nu0*kt)/(1 - nu*kt)
syms R kt nus
Eq_SG = R * (1 - nu0*kt) / (1 - nus*kt);
dSG_dR1  = double(subs(diff(Eq_SG, R),  [R,kt,nus],[R_sg1,kt_sg1,nu_spec]));
dSG_dkt1 = double(subs(diff(Eq_SG, kt), [R,kt,nus],[R_sg1,kt_sg1,nu_spec]));
dSG_dnu1 = double(subs(diff(Eq_SG, nus),[R,kt,nus],[R_sg1,kt_sg1,nu_spec]));
u_eps_sg1 = sqrt( (dSG_dR1*u_sg1_raw)^2 + (dSG_dkt1*s_kt_sg1_abs)^2 + (dSG_dnu1*s_nu)^2 );

dSG_dR2  = double(subs(diff(Eq_SG, R),  [R,kt,nus],[R_sg2,kt_sg2,nu_spec]));
dSG_dkt2 = double(subs(diff(Eq_SG, kt), [R,kt,nus],[R_sg2,kt_sg2,nu_spec]));
dSG_dnu2 = double(subs(diff(Eq_SG, nus),[R,kt,nus],[R_sg2,kt_sg2,nu_spec]));
u_eps_sg2 = sqrt( (dSG_dR2*u_sg2_raw)^2 + (dSG_dkt2*s_kt_sg2_abs)^2 + (dSG_dnu2*s_nu)^2 );

% Uniaxial stress & uncertainty: σ = E * ε_l
sigma_sg1 = E * eps_sg1;
sigma_sg2 = E * eps_sg2;
u_sigma_sg1 = sqrt( (eps_sg1*sE)^2 + (E*u_eps_sg1)^2 ); % dσ/dE = ε, dσ/dε = E
u_sigma_sg2 = sqrt( (eps_sg2*sE)^2 + (E*u_eps_sg2)^2 );

%% 6) ANGLE GAUGE (tee-rosette) — corrected strains & uncertainties
% Micro-Measurements form:
% ε_xx = ((1 - nu0*kt)*(εx - kt*εy)) / (1 - kt^2)
% ε_yy = ((1 - nu0*kt)*(εy - kt*εx)) / (1 - kt^2)
eps_ang_xx = ((1 - nu0*kt_ang) * (R_ang_x - kt_ang*R_ang_y)) / (1 - kt_ang^2);
eps_ang_yy = ((1 - nu0*kt_ang) * (R_ang_y - kt_ang*R_ang_x)) / (1 - kt_ang^2);

% Sensitivities
syms ex ey kt
Eq_xx = ((1 - nu0*kt) * (ex - kt*ey)) / (1 - kt^2);
Eq_yy = ((1 - nu0*kt) * (ey - kt*ex)) / (1 - kt^2);
dxx_dex = double(subs(diff(Eq_xx, ex), [ex,ey,kt],[R_ang_x,R_ang_y,kt_ang]));
dxx_dey = double(subs(diff(Eq_xx, ey), [ex,ey,kt],[R_ang_x,R_ang_y,kt_ang]));
dxx_dkt = double(subs(diff(Eq_xx, kt), [ex,ey,kt],[R_ang_x,R_ang_y,kt_ang]));
dyy_dex = double(subs(diff(Eq_yy, ex), [ex,ey,kt],[R_ang_x,R_ang_y,kt_ang]));
dyy_dey = double(subs(diff(Eq_yy, ey), [ex,ey,kt],[R_ang_x,R_ang_y,kt_ang]));
dyy_dkt = double(subs(diff(Eq_yy, kt), [ex,ey,kt],[R_ang_x,R_ang_y,kt_ang]));

% Uncertainties of corrected strains (RSS)
u_eps_ang_xx = sqrt( (dxx_dex*u_ang_xraw)^2 + (dxx_dey*u_ang_yraw)^2 + (dxx_dkt*s_kt_ang_abs)^2 );
u_eps_ang_yy = sqrt( (dyy_dex*u_ang_xraw)^2 + (dyy_dey*u_ang_yraw)^2 + (dyy_dkt*s_kt_ang_abs)^2 );

% Plane-stress (specimen Poisson nu_spec)
D = 1 - nu_spec^2;
sigma_ang_x = E/D * (eps_ang_xx + nu_spec*eps_ang_yy);
sigma_ang_y = E/D * (eps_ang_yy + nu_spec*eps_ang_xx);

% Stress sensitivities
dSx_dE   = (eps_ang_xx + nu_spec*eps_ang_yy)/D;
dSx_dexx = E/D;
dSx_deyy = E*nu_spec/D;
dSx_dnu  = E*( eps_ang_yy/D + (eps_ang_xx + nu_spec*eps_ang_yy)*(2*nu_spec)/D^2 );
dSy_dE   = (eps_ang_yy + nu_spec*eps_ang_xx)/D;
dSy_deyy = E/D;
dSy_dexx = E*nu_spec/D;
dSy_dnu  = E*( eps_ang_xx/D + (eps_ang_yy + nu_spec*eps_ang_xx)*(2*nu_spec)/D^2 );

u_sigma_ang_x = sqrt( (dSx_dE*sE)^2 + (dSx_dexx*u_eps_ang_xx)^2 + (dSx_deyy*u_eps_ang_yy)^2 + (dSx_dnu*s_nu)^2 );
u_sigma_ang_y = sqrt( (dSy_dE*sE)^2 + (dSy_dexx*u_eps_ang_xx)^2 + (dSy_deyy*u_eps_ang_yy)^2 + (dSy_dnu*s_nu)^2 );

%% 7) ROSETTE (0/45/90) — corrected ε_x, ε_y, γ_xy and uncertainties
% Normal components
eps_ros_xx = ((1 - nu0*kt_ros) * (R_ros_0  - kt_ros*R_ros_90)) / (1 - kt_ros^2);
eps_ros_yy = ((1 - nu0*kt_ros) * (R_ros_90 - kt_ros*R_ros_0 )) / (1 - kt_ros^2);

% Sensitivities for eps_ros_xx and eps_ros_yy
syms r0 r90 k
Eq_Rxx  = ((1 - nu0*k) * (r0  - k*r90)) / (1 - k^2);
Eq_Ryy  = ((1 - nu0*k) * (r90 - k*r0 )) / (1 - k^2);
drxx_dr0 = double(subs(diff(Eq_Rxx, r0),  [r0,r90,k],[R_ros_0,R_ros_90,kt_ros]));
drxx_dr9 = double(subs(diff(Eq_Rxx, r90), [r0,r90,k],[R_ros_0,R_ros_90,kt_ros]));
drxx_dk  = double(subs(diff(Eq_Rxx, k),   [r0,r90,k],[R_ros_0,R_ros_90,kt_ros]));
dryy_dr0 = double(subs(diff(Eq_Ryy, r0),  [r0,r90,k],[R_ros_0,R_ros_90,kt_ros]));
dryy_dr9 = double(subs(diff(Eq_Ryy, r90), [r0,r90,k],[R_ros_0,R_ros_90,kt_ros]));
dryy_dk  = double(subs(diff(Eq_Ryy, k),   [r0,r90,k],[R_ros_0,R_ros_90,kt_ros]));

u_eps_ros_xx = sqrt( (drxx_dr0*u_ros_0raw)^2 + (drxx_dr9*u_ros_90raw)^2 + (drxx_dk*s_kt_ros_abs)^2 );
u_eps_ros_yy = sqrt( (dryy_dr0*u_ros_0raw)^2 + (dryy_dr9*u_ros_90raw)^2 + (dryy_dk*s_kt_ros_abs)^2 );

% Shear: (kt_ros=0 → factor = 1)
gamma_raw   = 2*R_ros_45 - (R_ros_0 + R_ros_90);
u_gamma_raw = sqrt( (2*u_ros_45raw)^2 + (u_ros_0raw)^2 + (u_ros_90raw)^2 );

C_gamma = (1 - nu0*kt_ros) / (1 - kt_ros); % equals 1 for kt_ros=0
dC_dkt  = ( -(nu0)*(1 - kt_ros) + (1 - nu0*kt_ros) ) / (1 - kt_ros)^2;

gamma_ros_corr = C_gamma * gamma_raw;
u_gamma_ros    = sqrt( (C_gamma * u_gamma_raw)^2 + (dC_dkt * gamma_raw * s_kt_ros_abs)^2 );

% Plane-stress (specimen)
D = 1 - nu_spec^2;
sigma_ros_x = E/D * (eps_ros_xx + nu_spec*eps_ros_yy);
sigma_ros_y = E/D * (eps_ros_yy + nu_spec*eps_ros_xx);

% G = E/(1+nu_spec) 
tau_ros_xy = E/(2*(1 + nu_spec)) * gamma_ros_corr;

% Stress uncertainties (normals)
dSxR_dE   = (eps_ros_xx + nu_spec*eps_ros_yy)/D;
dSxR_dexx = E/D;
dSxR_deyy = E*nu_spec/D;
dSxR_dnu  = E*( eps_ros_yy/D + (eps_ros_xx + nu_spec*eps_ros_yy)*(2*nu_spec)/D^2 );
u_sigma_ros_x = sqrt( (dSxR_dE*sE)^2 + (dSxR_dexx*u_eps_ros_xx)^2 + (dSxR_deyy*u_eps_ros_yy)^2 + (dSxR_dnu*s_nu)^2 );

dSyR_dE   = (eps_ros_yy + nu_spec*eps_ros_xx)/D;
dSyR_deyy = E/D;
dSyR_dexx = E*nu_spec/D;
dSyR_dnu  = E*( eps_ros_xx/D + (eps_ros_yy + nu_spec*eps_ros_xx)*(2*nu_spec)/D^2 );
u_sigma_ros_y = sqrt( (dSyR_dE*sE)^2 + (dSyR_dexx*u_eps_ros_xx)^2 + (dSyR_deyy*u_eps_ros_yy)^2 + (dSyR_dnu*s_nu)^2 );

% Shear stress uncertainty -----
dT_dE = gamma_ros_corr/(2*(1 + nu_spec));
dT_dnu = -E * gamma_ros_corr / (2*(1 + nu_spec)^2);
dT_dg = E/(2*(1 + nu_spec));

u_tau_ros = sqrt( (dT_dE*sE)^2 + (dT_dnu*s_nu)^2 + (dT_dg*u_gamma_ros)^2 );

%% 8) RESULTS
fprintf('================ RESULTS ================\n\n');

% --- Single gauge ---
fprintf('SINGLE GAUGE\n');
fprintf(' ε_l : %8.3f ± %6.3f µε\n', 1e6*eps_sg1, 1e6*u_eps_sg1);
fprintf(' σ1 : %8.3f ± %6.3f MPa\n\n', sigma_sg1, u_sigma_sg1);
fprintf(' ε_2 : %8.3f ± %6.3f µε\n', 1e6*eps_sg2, 1e6*u_eps_sg2);
fprintf(' σ2 : %8.3f ± %6.3f MPa\n\n', sigma_sg2, u_sigma_sg2);

% --- Angle gauge ---
fprintf('ANGLE GAUGE (tee-rosette)\n');
fprintf(' ε_xx : %8.3f ± %6.3f µε\n', 1e6*eps_ang_xx, 1e6*u_eps_ang_xx);
fprintf(' ε_yy : %8.3f ± %6.3f µε\n', 1e6*eps_ang_yy, 1e6*u_eps_ang_yy);
fprintf(' σ_x : %8.3f ± %6.3f MPa\n', sigma_ang_x, u_sigma_ang_x);
fprintf(' σ_y : %8.3f ± %6.3f MPa\n\n', sigma_ang_y, u_sigma_ang_y);

% --- Rosette ---
fprintf('ROSETTE (0/45/90)\n');
fprintf(' ε_xx : %8.3f ± %6.3f µε\n', 1e6*eps_ros_xx, 1e6*u_eps_ros_xx);
fprintf(' ε_yy : %8.3f ± %6.3f µε\n', 1e6*eps_ros_yy, 1e6*u_eps_ros_yy);
fprintf(' γ_xy : %8.3f ± %6.3f µε\n', 1e6*gamma_ros_corr, 1e6*u_gamma_ros);
fprintf(' σ_x : %8.3f ± %6.3f MPa\n', sigma_ros_x, u_sigma_ros_x);
fprintf(' σ_y : %8.3f ± %6.3f MPa\n', sigma_ros_y, u_sigma_ros_y);
fprintf(' τ_xy : %8.3f ± %6.3f MPa\n', tau_ros_xy, u_tau_ros);
fprintf('=======================================================\n');

% Relative uncertainties
rel = @(val,unc) 100*unc/abs(val);
fprintf('\nRELATIVE UNCERTAINTIES (%%)\n');
fprintf(' Single: ε_l: %.3f%%, σ1: %.3f%%\n', rel(eps_sg1,u_eps_sg1), rel(sigma_sg1,u_sigma_sg1));
fprintf(' Single: ε_2: %.3f%%, σ2: %.3f%%\n', rel(eps_sg2,u_eps_sg2), rel(sigma_sg2,u_sigma_sg2));
fprintf(' Angle : ε_xx: %.3f%%, ε_yy: %.3f%%, σ_x: %.3f%%, σ_y: %.3f%%\n', ...
    rel(eps_ang_xx,u_eps_ang_xx), rel(eps_ang_yy,u_eps_ang_yy), rel(sigma_ang_x,u_sigma_ang_x), rel(sigma_ang_y,u_sigma_ang_y));
fprintf(' Roset : ε_xx: %.3f%%, ε_yy: %.3f%%, γ_xy: %.3f%%, σ_x: %.3f%%, σ_y: %.3f%%, τ_xy: %.3f%%\n', ...
    rel(eps_ros_xx,u_eps_ros_xx), rel(eps_ros_yy,u_eps_ros_yy), rel(gamma_ros_corr,u_gamma_ros), ...
    rel(sigma_ros_x,u_sigma_ros_x), rel(sigma_ros_y,u_sigma_ros_y), rel(tau_ros_xy,u_tau_ros));
