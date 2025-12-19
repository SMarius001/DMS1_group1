clc; clear; close all;

% Fonts used
set(0, 'DefaultAxesFontName', 'cmu');
set(0, 'DefaultTextFontName', 'cmu');
set(0, 'DefaultColorbarFontName', 'cmu');
set(0, 'DefaultLegendFontName', 'cmu');

% Font sizes used
set(0, 'DefaultAxesFontSize', 25);
set(0, 'DefaultTextFontSize', 25);
set(0, 'DefaultColorbarFontSize', 25);
set(0, 'DefaultLegendFontSize', 25);

%% ============================================================
% Load Excel file (DIC data)
%% ============================================================
DICdata = 'DatafromDIC.xlsx';
Forces  = readmatrix(DICdata, 'Range', 'H4:H11');

% Convert force from N to kN for plotting (y-axis)
Forces = Forces / 1000;

%% ============================================================
% Create / clean output folder for figures
%% ============================================================
saveFolder = fullfile(pwd, 'Figures_DOTDATA_DIC');
if exist(saveFolder, 'dir')
    delete(fullfile(saveFolder, '*'));   % delete everything inside
else
    mkdir(saveFolder);
end

%% ------------------------------------------------------------
% Hole loops (DIC data)  - exx, eyy, gxy (engineering shear)
%% ------------------------------------------------------------
Hole1_eyy = readmatrix(DICdata, 'Range', 'I4:I11');
Hole1_exx = readmatrix(DICdata, 'Range', 'J4:J11');
Hole1_gxy = 2 * readmatrix(DICdata, 'Range', 'K4:K11');  % γxy

Hole2_eyy = readmatrix(DICdata, 'Range', 'I14:I21');
Hole2_exx = readmatrix(DICdata, 'Range', 'J14:J21');
Hole2_gxy = 2 * readmatrix(DICdata, 'Range', 'K14:K21');

Hole3_eyy = readmatrix(DICdata, 'Range', 'I24:I31');
Hole3_exx = readmatrix(DICdata, 'Range', 'J24:J31');
Hole3_gxy = 2 * readmatrix(DICdata, 'Range', 'K24:K31');

% Make all DIC hole εyy positive
Hole1_eyy = abs(Hole1_eyy);
Hole2_eyy = abs(Hole2_eyy);
Hole3_eyy = abs(Hole3_eyy);

%% ------------------------------------------------------------
% Curve SG2 (DIC data)
%% ------------------------------------------------------------
Curve1_eyy_SG2 = readmatrix(DICdata, 'Range', 'B4:B11');
Curve1_exx_SG2 = readmatrix(DICdata, 'Range', 'C4:C11');
Curve1_gxy_SG2 = 2 * readmatrix(DICdata, 'Range', 'D4:D11');

Curve2_eyy_SG2 = readmatrix(DICdata, 'Range', 'B14:B21');
Curve2_exx_SG2 = readmatrix(DICdata, 'Range', 'C14:C21');
Curve2_gxy_SG2 = 2 * readmatrix(DICdata, 'Range', 'D14:D21');

Curve3_eyy_SG2 = readmatrix(DICdata, 'Range', 'B24:B31');
Curve3_exx_SG2 = readmatrix(DICdata, 'Range', 'C24:C31');
Curve3_gxy_SG2 = 2 * readmatrix(DICdata, 'Range', 'D24:D31');

%% ------------------------------------------------------------
% Curve SG1 (DIC data)
%% ------------------------------------------------------------
Curve1_eyy_SG1 = readmatrix(DICdata, 'Range', 'P4:P11');
Curve1_exx_SG1 = readmatrix(DICdata, 'Range', 'Q4:Q11');
Curve1_gxy_SG1 = 2 * readmatrix(DICdata, 'Range', 'R4:R11');

Curve2_eyy_SG1 = readmatrix(DICdata, 'Range', 'P14:P21');
Curve2_exx_SG1 = readmatrix(DICdata, 'Range', 'Q14:Q21');  % FIXED (was Q14:P21)
Curve2_gxy_SG1 = 2 * readmatrix(DICdata, 'Range', 'R14:R21');

Curve3_eyy_SG1 = readmatrix(DICdata, 'Range', 'P24:P31');
Curve3_exx_SG1 = readmatrix(DICdata, 'Range', 'Q24:Q31');
Curve3_gxy_SG1 = 2 * readmatrix(DICdata, 'Range', 'R24:R31');

%% ============================================================
% Load strain gauge comparison data
%% ============================================================
SGdata  = 'DataComparisonStrainGaugesandForce.xlsx';
sgSheet = 'Comparison and Org';

% Hole strain gauge (εxx, εyy, ε45)
Hole_eyy_SG = readGaugeStrain(SGdata, sgSheet, 'AO');
Hole_exx_SG = readGaugeStrain(SGdata, sgSheet, 'AN');
Hole_e45_SG = readGaugeStrain(SGdata, sgSheet, 'AM');

% γxy (engineering shear) = 2ε45 − εxx − εyy
Hole_gxy_SG = 2*Hole_e45_SG - Hole_exx_SG - Hole_eyy_SG;

% Curve SG2 and SG1 εyy equivalent
Curve_SG2_eyy = readGaugeStrain(SGdata, sgSheet, 'AR');
Curve_SG1_eyy = readGaugeStrain(SGdata, sgSheet, 'AS');

% Corresponding force levels for the strain gauge data [N]
Forces_SG = (500:500:4000).';

% Convert force from N to kN for plotting (y-axis)
Forces_SG = Forces_SG / 1000;

%% ============================================================
% Material properties
%% ============================================================
E  = 205e3;   % MPa
nu = 0.253;
D  = (E/(1 - nu^2)) * [ 1   nu   0;
                        nu  1    0;
                        0   0   (1 - nu)/2 ];

%% ============================================================
% Compute mean strains (engineering γxy)
%% ============================================================
Hole_EYY_mean = mean([Hole1_eyy, Hole2_eyy, Hole3_eyy], 2);
Hole_EXX_mean = mean([Hole1_exx, Hole2_exx, Hole3_exx], 2);
Hole_GXY_mean = mean([Hole1_gxy, Hole2_gxy, Hole3_gxy], 2);

SG2_EYY_mean  = mean([Curve1_eyy_SG2, Curve2_eyy_SG2, Curve3_eyy_SG2], 2);
SG2_EXX_mean  = mean([Curve1_exx_SG2, Curve2_exx_SG2, Curve3_exx_SG2], 2);
SG2_GXY_mean  = mean([Curve1_gxy_SG2, Curve2_gxy_SG2, Curve3_gxy_SG2], 2);

SG1_EYY_mean  = mean([Curve1_eyy_SG1, Curve2_eyy_SG1, Curve3_eyy_SG1], 2);
SG1_EXX_mean  = mean([Curve1_exx_SG1, Curve2_exx_SG1, Curve3_exx_SG1], 2);
SG1_GXY_mean  = mean([Curve1_gxy_SG1, Curve2_gxy_SG1, Curve3_gxy_SG1], 2);

%% ============================================================
% Strain & stress naming conventions
%% ============================================================
strain_names = {'\epsilon_{xx}','\epsilon_{yy}','\gamma_{xy}'};
strain_tags  = {'eps_xx','eps_yy','gamma_xy'};

stress_names = {'\sigma_{xx}','\sigma_{yy}','\tau_{xy}'};
stress_tags  = {'sig_xx','sig_yy','tau_xy'};

%% ============================================================
% ------------------------- HOLE STRAIN & STRESS -------------------------
%% ============================================================
Hole_data = {
    Hole1_exx,     Hole1_eyy,     Hole1_gxy;
    Hole2_exx,     Hole2_eyy,     Hole2_gxy;
    Hole3_exx,     Hole3_eyy,     Hole3_gxy;
    Hole_EXX_mean, Hole_EYY_mean, Hole_GXY_mean
};

Hole_labels = {'Hole loop 1','Hole loop 2','Hole loop 3','Mean'};

% ----- Strain plots for holes (εxx, εyy, γxy) -----
for comp = 1:3
    baseName = sprintf('Hole_strain_%s_vs_Force', strain_tags{comp});
    fig = figure('Name', baseName, 'NumberTitle', 'off');
    hold on; grid on; box on;

    for dataset = 1:4
        plot(Hole_data{dataset, comp}, Forces, '-o', 'LineWidth', 1.6);
    end

    switch comp
        case 1
            plot(Hole_exx_SG, Forces_SG, 's--', 'LineWidth', 1.6);
        case 2
            plot(Hole_eyy_SG, Forces_SG, 's--', 'LineWidth', 1.6);
        case 3
            plot(Hole_gxy_SG, Forces_SG, 's--', 'LineWidth', 1.6); % γxy
    end

    xlabel([strain_names{comp}]);
    ylabel('Force [kN]');
    legend([Hole_labels, {'Strain gauge'}], 'Location','best');
    hold off;

    exportgraphics(fig, fullfile(saveFolder,[baseName '.png']), 'Resolution',300);
end

% ----- Compute stresses (DIC + SG) for holes -----
Hole_stress_all    = cell(1,4);
Hole_vm_strain_all = cell(1,4);
Hole_vm_stress_all = cell(1,4);

for j = 1:4
    [s, vm_e, vm_s] = computeStressVM( ...
        Hole_data{j,1}, Hole_data{j,2}, Hole_data{j,3}, D);
    Hole_stress_all{j}    = s;
    Hole_vm_strain_all{j} = vm_e;
    Hole_vm_stress_all{j} = vm_s;
end

[Hole_stress_SG, Hole_vm_strain_SG, Hole_vm_stress_SG] = computeStressVM( ...
    Hole_exx_SG, Hole_eyy_SG, Hole_gxy_SG, D);

% ----- Stress plots for holes (σxx, σyy, τxy) -----
for comp = 1:3
    baseName = sprintf('Hole_stress_%s_vs_Force', stress_tags{comp});
    fig = figure('Name', baseName, 'NumberTitle', 'off');
    hold on; grid on; box on;

    for dataset = 1:4
        plot(Hole_stress_all{dataset}(:, comp), Forces, '-o', 'LineWidth',1.6);
    end

    plot(Hole_stress_SG(:, comp), Forces_SG, 's--', 'LineWidth', 1.6);

    xlabel([stress_names{comp} ' [MPa]']);
    ylabel('Force [kN]');
    legend([Hole_labels, {'Strain gauge'}],'Location','best');
    hold off;

    exportgraphics(fig, fullfile(saveFolder,[baseName '.png']), 'Resolution',300);
end

% ----- Equivalent strain (holes) -----
baseName = 'Hole_eq_strain_vm_vs_Force';
fig = figure('Name', baseName, 'NumberTitle', 'off');
hold on; grid on; box on;
for dataset = 1:4
    plot(Hole_vm_strain_all{dataset}, Forces, '-o', 'LineWidth', 1.6);
end
plot(Hole_vm_strain_SG, Forces_SG, 's--', 'LineWidth', 1.6);
xlabel('\epsilon_{vm}');
ylabel('Force [kN]');
legend([Hole_labels, {'Strain gauge'}], 'Location', 'best');
hold off;
exportgraphics(fig, fullfile(saveFolder, [baseName '.png']), 'Resolution', 300);

% ----- Equivalent stress (holes) -----
baseName = 'Hole_eq_stress_vm_vs_Force';
fig = figure('Name', baseName, 'NumberTitle', 'off');
hold on; grid on; box on;
for dataset = 1:4
    plot(Hole_vm_stress_all{dataset}, Forces, '-o', 'LineWidth', 1.6);
end
plot(Hole_vm_stress_SG, Forces_SG, 's--', 'LineWidth', 1.6);
xlabel('\sigma_{eq} [MPa]');
ylabel('Force [kN]');
legend([Hole_labels, {'Strain gauge'}], 'Location', 'best');
hold off;
exportgraphics(fig, fullfile(saveFolder, [baseName '.png']), 'Resolution', 300);

%% ============================================================
% ------------------------- CURVE SG2 (DIC at SG2) -------------------------
%% ============================================================
Curve_data = {
    Curve1_exx_SG2, Curve1_eyy_SG2, Curve1_gxy_SG2;
    Curve2_exx_SG2, Curve2_eyy_SG2, Curve2_gxy_SG2;
    Curve3_exx_SG2, Curve3_eyy_SG2, Curve3_gxy_SG2;
    SG2_EXX_mean,   SG2_EYY_mean,   SG2_GXY_mean
};

Curve_labels = {'Curve 1 SG2','Curve 2 SG2','Curve 3 SG2','Mean SG2'};

% ----- Strain plots for Curve SG2 -----
for comp = 1:3
    baseName = sprintf('CurveSG2_strain_%s_vs_Force', strain_tags{comp});
    fig = figure('Name', baseName, 'NumberTitle', 'off');
    hold on; grid on; box on;

    for dataset = 1:4
        plot(Curve_data{dataset, comp}, Forces, '-o', 'LineWidth', 1.6);
    end

    if comp == 2
        plot(Curve_SG2_eyy, Forces_SG, 's--', 'LineWidth', 1.6);
        legend_entries = [Curve_labels, {'Strain gauge'}];
    else
        legend_entries = Curve_labels;
    end

    xlabel([strain_names{comp}]);
    ylabel('Force [kN]');
    legend(legend_entries, 'Location', 'best');
    hold off;

    exportgraphics(fig, fullfile(saveFolder, [baseName '.png']), 'Resolution', 300);
end

% ----- Stresses & vm for Curve SG2 (DIC) -----
Curve_stress_all     = cell(1,4);
Curve_vm_strain_all  = cell(1,4);
Curve_vm_stress_all  = cell(1,4);

for j = 1:4
    [s, vm_e, vm_s] = computeStressVM( ...
        Curve_data{j,1}, Curve_data{j,2}, Curve_data{j,3}, D);
    Curve_stress_all{j}    = s;
    Curve_vm_strain_all{j} = vm_e;
    Curve_vm_stress_all{j} = vm_s;
end

% Strain gauge stress for SG2 (assume exx=0, γxy=0)
zero_SG2 = zeros(size(Curve_SG2_eyy));
[Curve_SG2_stress_SG, Curve_SG2_vm_strain_SG, Curve_SG2_vm_stress_SG] = ...
    computeStressVM(zero_SG2, Curve_SG2_eyy, zero_SG2, D);

% ----- Stress plots for Curve SG2 -----
for comp = 1:3
    baseName = sprintf('CurveSG2_stress_%s_vs_Force', stress_tags{comp});
    fig = figure('Name', baseName, 'NumberTitle', 'off');
    hold on; grid on; box on;

    for dataset = 1:4
        plot(Curve_stress_all{dataset}(:, comp), Forces, '-o', 'LineWidth', 1.6);
    end

    if comp == 2
        plot(Curve_SG2_stress_SG(:, 2), Forces_SG, 's--', 'LineWidth', 1.6);
        legend_entries = [Curve_labels, {'Strain gauge'}];
    else
        legend_entries = Curve_labels;
    end

    xlabel([stress_names{comp} ' [MPa]']);
    ylabel('Force [kN]');
    legend(legend_entries, 'Location', 'best');
    hold off;

    exportgraphics(fig, fullfile(saveFolder, [baseName '.png']), 'Resolution', 300);
end

% ----- Equivalent strain (Curve SG2, DIC only) -----
baseName = 'CurveSG2_eq_strain_vm_vs_Force';
fig = figure('Name', baseName, 'NumberTitle', 'off');
hold on; grid on; box on;
for dataset = 1:4
    plot(Curve_vm_strain_all{dataset}, Forces, '-o', 'LineWidth', 1.6);
end
xlabel('\epsilon_{eq}');
ylabel('Force [kN]');
legend(Curve_labels, 'Location', 'best');
hold off;
exportgraphics(fig, fullfile(saveFolder, [baseName '.png']), 'Resolution', 300);

% ----- Equivalent stress (Curve SG2, DIC only) -----
baseName = 'CurveSG2_eq_stress_vm_vs_Force';
fig = figure('Name', baseName, 'NumberTitle', 'off');
hold on; grid on; box on;
for dataset = 1:4
    plot(Curve_vm_stress_all{dataset}, Forces, '-o', 'LineWidth', 1.6);
end
xlabel('\sigma_{eq} [MPa]');
ylabel('Force [kN]');
legend(Curve_labels, 'Location', 'best');
hold off;
exportgraphics(fig, fullfile(saveFolder, [baseName '.png']), 'Resolution', 300);

%% ============================================================
% ------------------------- CURVE SG1 (DIC at SG1) -------------------------
% ============================================================
Curve2_data = {
    Curve1_exx_SG1, Curve1_eyy_SG1, Curve1_gxy_SG1;
    Curve2_exx_SG1, Curve2_eyy_SG1, Curve2_gxy_SG1;
    Curve3_exx_SG1, Curve3_eyy_SG1, Curve3_gxy_SG1;
    SG1_EXX_mean,   SG1_EYY_mean,   SG1_GXY_mean
};

Curve2_labels = {'Curve 1 SG1','Curve 2 SG1','Curve 3 SG1','Mean SG1'};

% ----- Strain plots for Curve SG1 -----
for comp = 1:3
    baseName = sprintf('CurveSG1_strain_%s_vs_Force', strain_tags{comp});
    fig = figure('Name', baseName, 'NumberTitle', 'off');
    hold on; grid on; box on;

    for dataset = 1:4
        plot(Curve2_data{dataset, comp}, Forces, '-o', 'LineWidth', 1.6);
    end

    if comp == 2
        plot(Curve_SG1_eyy, Forces_SG, 's--', 'LineWidth', 1.6);
        legend_entries = [Curve2_labels, {'Strain gauge'}];
    else
        legend_entries = Curve2_labels;
    end

    xlabel([strain_names{comp}]);
    ylabel('Force [kN]');
    legend(legend_entries, 'Location', 'best');
    hold off;

    exportgraphics(fig, fullfile(saveFolder, [baseName '.png']), 'Resolution', 300);
end

% ----- Stresses & vm for Curve SG1 (DIC) -----
Curve2_stress_all     = cell(1,4);
Curve2_vm_strain_all  = cell(1,4);
Curve2_vm_stress_all  = cell(1,4);

for j = 1:4
    [s, vm_e, vm_s] = computeStressVM( ...
        Curve2_data{j,1}, Curve2_data{j,2}, Curve2_data{j,3}, D);
    Curve2_stress_all{j}    = s;
    Curve2_vm_strain_all{j} = vm_e;
    Curve2_vm_stress_all{j} = vm_s;
end

% Strain gauge stress for SG1 (assume exx=0, γxy=0)
zero_SG1 = zeros(size(Curve_SG1_eyy));
[Curve_SG1_stress_SG, Curve_SG1_vm_strain_SG, Curve_SG1_vm_stress_SG] = ...
    computeStressVM(zero_SG1, Curve_SG1_eyy, zero_SG1, D);

% ----- Stress plots for Curve SG1 -----
for comp = 1:3
    baseName = sprintf('CurveSG1_stress_%s_vs_Force', stress_tags{comp});
    fig = figure('Name', baseName, 'NumberTitle', 'off');
    hold on; grid on; box on;

    for dataset = 1:4
        plot(Curve2_stress_all{dataset}(:, comp), Forces, '-o', 'LineWidth', 1.6);
    end

    if comp == 2
        plot(Curve_SG1_stress_SG(:, 2), Forces_SG, 's--', 'LineWidth', 1.6);
        legend_entries = [Curve2_labels, {'Strain gauge'}];
    else
        legend_entries = Curve2_labels;
    end

    xlabel([stress_names{comp} ' [MPa]']);
    ylabel('Force [kN]');
    legend(legend_entries, 'Location', 'best');
    hold off;

    exportgraphics(fig, fullfile(saveFolder, [baseName '.png']), 'Resolution', 300);
end

% ----- Equivalent strain (Curve SG1, DIC only) -----
baseName = 'CurveSG1_eq_strain_vm_vs_Force';
fig = figure('Name', baseName, 'NumberTitle', 'off');
hold on; grid on; box on;
for dataset = 1:4
    plot(Curve2_vm_strain_all{dataset}, Forces, '-o', 'LineWidth', 1.6);
end
xlabel('\epsilon_{eq}');
ylabel('Force [kN]');
legend(Curve2_labels, 'Location', 'best');
hold off;
exportgraphics(fig, fullfile(saveFolder, [baseName '.png']), 'Resolution', 300);

% ----- Equivalent stress (Curve SG1, DIC only) -----
baseName = 'CurveSG1_eq_stress_vm_vs_Force';
fig = figure('Name', baseName, 'NumberTitle', 'off');
hold on; grid on; box on;
for dataset = 1:4
    plot(Curve2_vm_stress_all{dataset}, Forces, '-o', 'LineWidth', 1.6);
end
xlabel('\sigma_{eq} [MPa]');
ylabel('Force [kN]');
legend(Curve2_labels, 'Location', 'best');
hold off;
exportgraphics(fig, fullfile(saveFolder, [baseName '.png']), 'Resolution', 300);

%% ============================================================
% ------------------------- FUNCTIONS -------------------------
%% ============================================================
function [stress, vm_strain, vm_stress] = computeStressVM(exx, eyy, gxy, D)
    % exx, eyy, gxy are ENGINEERING strains (γ_xy).
    % Strain vector for plane stress: [εxx, εyy, γxy]
    strain_matrix = [exx, eyy, gxy];

    % Stress using constitutive matrix D: [σxx, σyy, τxy]
    stress = (D * strain_matrix')';

    % Tensor shear strain E_xy = γxy / 2
    E_xy = gxy / 2;

    % von Mises equivalent strain (based on tensor components)
    vm_strain = sqrt(exx.^2 - exx .* eyy + eyy.^2 + 3 * E_xy.^2);

    % von Mises equivalent stress
    vm_stress = sqrt( ...
        stress(:,1).^2 - stress(:,1).*stress(:,2) + stress(:,2).^2 + ...
        3 * stress(:,3).^2 );
end

function strainVals = readGaugeStrain(filename, sheetName, colLetter)
    % Read ONLY the eight specified cells in the given column:
    % rows 141, 149, 157, 165, 173, 181, 189, 197
    rows = [141 149 157 165 173 181 189 197];

    n    = numel(rows);
    vals = nan(n,1);

    for k = 1:n
        cellRange = sprintf('%s%d:%s%d', colLetter, rows(k), colLetter, rows(k));
        tmp = readmatrix(filename, 'Sheet', sheetName, 'Range', cellRange);

        if isempty(tmp)
            vals(k) = NaN;
        else
            vals(k) = tmp(1);
        end
    end

    % Convert microstrain -> regular strain
    strainVals = vals * 1e-6;
end
