clc; clear; close all;

%% ============================================================
% Global font settings
%% ============================================================
% Fonts used
set(0, 'DefaultAxesFontName',    'cmu');
set(0, 'DefaultTextFontName',    'cmu');
set(0, 'DefaultColorbarFontName','cmu');
set(0, 'DefaultLegendFontName',  'cmu');

% Font sizes used
set(0, 'DefaultAxesFontSize',    18);
set(0, 'DefaultTextFontSize',    18);
set(0, 'DefaultColorbarFontSize',18);
set(0, 'DefaultLegendFontSize',  18);

%% ============================================================
% Create / clean output folder for all figures (vector graphics)
%% ============================================================
saveFolder = fullfile(pwd, 'Figures_finalpull_dotdata');

if exist(saveFolder, 'dir')
    delete(fullfile(saveFolder, '*'));   % delete everything inside
else
    mkdir(saveFolder);
end

%% ============================================================
% 1) Load Ncorr results (.mat file with DIC data)
%% ============================================================

ncorrFile = 'Finalpull_curve_Ncorrfile.mat';
raw = load(ncorrFile);

dicStruct = [];

if isfield(raw, 'handles_ncorr') && isfield(raw.handles_ncorr, 'data_dic')
    dicStruct = raw.handles_ncorr.data_dic;
end
if isempty(dicStruct) && isfield(raw, 'data_dic')
    dicStruct = raw.data_dic;
end
if isempty(dicStruct) && isfield(raw, 'data_dic_save')
    dicStruct = raw.data_dic_save;
end
if isempty(dicStruct)
    fns = fieldnames(raw);
    for i = 1:numel(fns)
        v = raw.(fns{i});
        if isstruct(v)
            if isfield(v, 'data_dic')
                dicStruct = v.data_dic;
                break;
            elseif isfield(v, 'strains')
                dicStruct = v;
                break;
            end
        end
    end
end

if isempty(dicStruct)
    error('Could not find a DIC structure in file "%s".', ncorrFile);
end

%% ============================================================
% 2) Choose coordinate of interest
%% ============================================================

x_pos = 245;   
y_pos = 176;

%% ============================================================
% 3) Extract strain at that point over all frames
%% ============================================================

useEulerian = false;
nFrames     = numel(dicStruct.strains);

Exx = nan(nFrames,1);
Eyy = nan(nFrames,1);
Exy = nan(nFrames,1);   % TO BE CONVERTED TO ENGINEERING FORM

for k = 1:nFrames
    if ~useEulerian
        exxField = dicStruct.strains(k).plot_exx_ref_formatted;
        eyyField = dicStruct.strains(k).plot_eyy_ref_formatted;
        exyField = dicStruct.strains(k).plot_exy_ref_formatted;
    else
        exxField = dicStruct.strains(k).plot_exx_cur_formatted;
        eyyField = dicStruct.strains(k).plot_eyy_cur_formatted;
        exyField = dicStruct.strains(k).plot_exy_cur_formatted;
    end

    if isempty(exxField)
        continue;
    end

    [ny, nx] = size(exxField);
    if y_pos <= ny && x_pos <= nx
        Exx(k) = exxField(y_pos, x_pos);
        Eyy(k) = eyyField(y_pos, x_pos);
        Exy(k) = exyField(y_pos, x_pos);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATCH: Convert DIC tensor shear strain → engineering shear strain
%        γ_xy = 2 * E_xy
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Exy = 2 * Exy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% ============================================================
% 4) Build force vector (0 → 16.86 kN)
%% ============================================================

F_max_kN = 16.86;
F_max_N  = F_max_kN * 1e3;

Force = linspace(0, F_max_N, nFrames).';

%% ============================================================
% 4b) Load comparison strain gauge data
%% ============================================================

compFile = 'DataComparisonStrainGaugesandForce.xlsx';

gaugeMicro = readmatrix(compFile, ...
    'Sheet', 'PlasticTestCurved', ...
    'Range', 'J997:J2300');

gaugeEyy = gaugeMicro * 1e-6;

ForceGauge = linspace(0, F_max_N, numel(gaugeEyy)).';

%% ============================================================
% 5) Plot strain vs force (NO TITLES)
%% ============================================================

% ---- Exx vs Force ----
fig = figure('Name','Exx vs Force');
plot(Force, Exx, '-o');
grid on;
xlabel('Force [N]');
ylabel('$\varepsilon_{xx}$','Interpreter','latex');
exportgraphics(fig, fullfile(saveFolder, 'Exx_vs_Force.pdf'), 'ContentType', 'vector');

% ---- Eyy vs Force ----
fig = figure('Name','Eyy vs Force');
plot(Force, Eyy, '-o');
hold on;
plot(ForceGauge, gaugeEyy, '-x');
hold off;
grid on;
xlabel('Force [N]');
ylabel('$\varepsilon_{yy}$','Interpreter','latex');
leg = legend({'DIC','Strain gauge'}, 'Location', 'best');
leg.FontSize = 16;
leg.ItemTokenSize = [30, 18];
exportgraphics(fig, fullfile(saveFolder, 'Eyy_vs_Force.pdf'), 'ContentType', 'vector');

% ---- Exy vs Force (ENGINEERING) ----
fig = figure('Name','Exy vs Force');
plot(Force, Exy, '-o');
grid on;
xlabel('Force [N]');
ylabel('$\gamma_{xy}$','Interpreter','latex');   % now engineering shear!
exportgraphics(fig, fullfile(saveFolder, 'Exy_vs_Force.pdf'), 'ContentType', 'vector');

%% ============================================================
% 6) Contour plots of strain fields
%% ============================================================

frameStep = 20;
framesToPlot = frameStep:frameStep:nFrames;

for idx = 1:numel(framesToPlot)
    k = framesToPlot(idx);

    if ~useEulerian
        exxField = dicStruct.strains(k).plot_exx_ref_formatted;
        eyyField = dicStruct.strains(k).plot_eyy_ref_formatted;
        exyField = dicStruct.strains(k).plot_exy_ref_formatted;
    else
        exxField = dicStruct.strains(k).plot_exx_cur_formatted;
        eyyField = dicStruct.strains(k).plot_eyy_cur_formatted;
        exyField = dicStruct.strains(k).plot_exy_cur_formatted;
    end

    if isempty(exxField)
        continue;
    end

    exxPlot = exxField;
    eyyPlot = eyyField;
    exyPlot = exyField;

    exxPlot(exxPlot == 0) = NaN;
    eyyPlot(eyyPlot == 0) = NaN;
    exyPlot(exyPlot == 0) = NaN;

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PATCH: Convert all DIC shear field values → engineering γxy
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    exyPlot = 2 * exyPlot;
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Determine cropping region
    validMask = ~(isnan(exxPlot) & isnan(eyyPlot) & isnan(exyPlot));
    if any(validMask(:))
        [rowIdx, colIdx] = find(validMask);
        rmin = min(rowIdx); rmax = max(rowIdx);
        cmin = min(colIdx); cmax = max(colIdx);
        crop = @(M) M(rmin:rmax, cmin:cmax);
    else
        crop = @(M) M;
        rmin = 1; rmax = size(exxPlot,1);
        cmin = 1; cmax = size(exxPlot,2);
    end

    if y_pos >= rmin && y_pos <= rmax && x_pos >= cmin && x_pos <= cmax
        dotRow = y_pos - rmin + 1;
        dotCol = x_pos - cmin + 1;
        showDot = true;
    else
        showDot = false;
    end

    dataCells   = {exxPlot, eyyPlot, exyPlot};
    labelStrain = {'\epsilon_{xx}','\epsilon_{yy}','\gamma_{xy}'}; % now engineering

    Fk_N  = Force(k);
    Fk_kN = Fk_N / 1e3;

    fig = figure('Name', sprintf('Strain fields at F = %.2f kN', Fk_kN));
    tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

    for j = 1:3
        ax = nexttile(tl);
        Z = crop(dataCells{j});

        contourf(ax, Z, 15, 'LineStyle', 'none');
        axis(ax, 'equal', 'tight');
        ax.XTick = [];
        ax.YTick = [];

        if showDot
            hold(ax, 'on');
            plot(ax, dotCol, dotRow, 'k.', 'MarkerSize', 20);
            hold(ax, 'off');
        end

        cb = colorbar(ax);
        cb.Label.Interpreter = 'latex';
        cb.Label.String = sprintf('$%s$ [-]', labelStrain{j});
        cb.Label.FontSize = 22;
        colormap(ax, jet);
    end

    set(fig, 'PaperPositionMode', 'auto');

    fileNameContour = sprintf('StrainFields_F_%.2fkN.pdf', Fk_kN);
    exportgraphics(fig, fullfile(saveFolder, fileNameContour), 'ContentType', 'vector');
end
