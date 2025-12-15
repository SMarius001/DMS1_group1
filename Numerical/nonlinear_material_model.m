%% SCRIPT TO CALCULATE TRUE STRESS-STRAIN CURVE FOR S355-J2
clear; clc; close all;

% Loading the experimental stress and strain
% sig_eps_eng = load('material_stress_strain_0.csv'); % For 0 deg
sig_eps_eng = load('material_stress_strain_90.csv'); % For 90 deg
sig_eps_eng(1:5,:) = []; % Removes elastic region
eps_yield_eng = sig_eps_eng(1,2) / 205e3;
sig_eps_eng(:,1) = sig_eps_eng(:,1) - (sig_eps_eng(1,1) - eps_yield_eng);
eps_eng = sig_eps_eng(:,1);
sig_eng = sig_eps_eng(:,2);

% Calculating the true stress and true strain
eps_true = log(1 + eps_eng);
sig_true = sig_eng .* exp(eps_true);
sig_eps_true = [eps_true sig_true];

% Calculate the plastic strain
eps_eng_elastic = sig_eng ./ 205e3;
eps_eng_plastic = eps_eng - eps_eng_elastic;
eps_true_elastic = sig_true ./ 205e3;
eps_true_plastic = eps_true - eps_true_elastic;

% Calculate the tangent modulus for a bi-linear hardening model
eps_yield_true = eps_true(1);
sig_yield_true = sig_true(1);
% eps_uts_true = eps_true(21); % For 0 deg
% sig_uts_true = sig_true(21); % For 0 deg
eps_uts_true = eps_true(18); % For 90 deg
sig_uts_true = sig_true(18); % For 90 deg
tan_modulus = (sig_uts_true - sig_yield_true) / (eps_uts_true - eps_yield_true);

% Define region to plot
% sig_eps_eng_plot = [eps_eng_plastic(1:21) sig_eng(1:21)]; % For 0 deg
% sig_eps_true_plot = [eps_true_plastic(1:21) sig_true(1:21)]; % For 0 deg
sig_eps_eng_plot = [eps_eng_plastic(1:18) sig_eng(1:18)]; % For 90 deg
sig_eps_true_plot = [eps_true_plastic(1:18) sig_true(1:18)]; % For 90 deg

%% Plot curves
figure
hold on
% plot(sig_eps_eng_plot(:,1),sig_eps_eng_plot(:,2),'DisplayName','Engineering')
plot(sig_eps_true_plot(:,1),sig_eps_true_plot(:,2),'DisplayName','True')
% plot([eps_true_plastic(1); eps_true_plastic(21)],[sig_true(1); sig_true(21)],'DisplayName','Bi-linear')
% plot([eps_true_plastic(1); eps_true_plastic(18)],[sig_true(1); sig_true(18)],'DisplayName','Bi-linear')
xlabel('Plastic Strain')
ylabel('Stress [MPa]')
% legend('Location','southeast')
% title('Plastic Stress-Strain Curve')
box on
grid on
axis padded
hold off
% set(gca, 'FontName', 'Helvetica')
set(gca, 'FontName', 'CMU Serif')
set(gca, 'FontSize', 10)
