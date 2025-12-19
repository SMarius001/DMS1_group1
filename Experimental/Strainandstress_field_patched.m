clc; clear; close all

% Fonts used
set(0, 'DefaultAxesFontName', 'cmu');
set(0, 'DefaultTextFontName', 'cmu');
set(0, 'DefaultColorbarFontName', 'cmu');
set(0, 'DefaultLegendFontName', 'cmu');

% Font sizes used
set(0, 'DefaultAxesFontSize', 18);
set(0, 'DefaultTextFontSize', 18);
set(0, 'DefaultColorbarFontSize', 18);
set(0, 'DefaultLegendFontSize', 18);

%% -------------------------------------------------------------
%  Create / clean output folder
%% -------------------------------------------------------------
saveFolder = fullfile(pwd, 'Figures_contourplots_DIC');

if exist(saveFolder, 'dir')
    delete(fullfile(saveFolder, '*'));   % delete everything inside
else
    mkdir(saveFolder);
end

%% -------------------------------------------------------------
%  Data Definitions
%% -------------------------------------------------------------
files = {
    'Hole1_data.xlsx'
    'Hole2_data.xlsx'
    'Hole3_data.xlsx'
    'Curveloop1_data.xlsx'
    'Curveloop2_data.xlsx'
    'Curveloop3_data.xlsx'
};

names = {
    'Hole loop 1'
    'Hole loop 2'
    'Hole loop 3'
    'Curve loop 1'
    'Curve loop 2'
    'Curve loop 3'
};

components = {'eyy', 'exx', 'exy'};

%% -------------------------------------------------------------
%  Material Parameters
%% -------------------------------------------------------------
E      = 205e9;           
E_MPa  = E / 1e6;         
nu     = 0.253;

% Shear modulus
G_MPa  = E_MPa / (2*(1 + nu));

%% -------------------------------------------------------------
%  Loop through datasets
%% -------------------------------------------------------------
for fi = 1:length(files)

    fileName    = files{fi};
    datasetName = names{fi};

    exx = []; eyy = []; exy = [];

    %% ---------------------------------------------------------
    % Load strain matrices
    %% ---------------------------------------------------------
    for ci = 1:length(components)

        comp      = components{ci};
        sheetName = sprintf('%s_4000', comp);

        strainMatrix = table2array(readtable(fileName, 'Sheet', sheetName));

        % Fix upside-down hole data (Hole 1–3)
        if fi <= 3
            strainMatrix = flipud(strainMatrix);
        end

        % Mask zeros as NaN
        strainMatrix(strainMatrix == 0) = NaN;

        switch comp
            case 'exx', exx = strainMatrix;
            case 'eyy', eyy = strainMatrix;
            case 'exy', exy = strainMatrix;
        end
    end

    %% ---------------------------------------------------------
    %  Compute 2D von Mises Stress
    %% ---------------------------------------------------------
    C       = E_MPa/(1 - nu^2);
    sigma_x = C * (exx + nu * eyy);
    sigma_y = C * (eyy + nu * exx);

    % -------------------------------------------------------------
    %  *** PATCHED SECTION: Correct shear stress calculation ***
    %  exy from Ncorr = tensor shear strain E_xy
    %  engineering shear strain = gamma_xy = 2*E_xy
    %  tau_xy = G * gamma_xy = 2*G*E_xy
    % -------------------------------------------------------------
    tau_xy  = 2 * G_MPa .* exy;
    % -------------------------------------------------------------

    % von Mises stress
    sigma_vm = sqrt( sigma_x.^2 + sigma_y.^2 - sigma_x .* sigma_y + 3 * tau_xy.^2 );

    maxVM = max(sigma_vm, [], 'all', 'omitnan');
    fprintf('Max von Mises stress in %s: %.3f MPa\n', datasetName, maxVM);

    %% ---------------------------------------------------------
    %  Determine cropping window
    %% ---------------------------------------------------------
    validMask = ~isnan(sigma_vm);
    if any(validMask(:))
        [rowIdx, colIdx] = find(validMask);
        rmin = min(rowIdx); rmax = max(rowIdx);
        cmin = min(colIdx); cmax = max(colIdx);

        crop = @(M) M(rmin:rmax, cmin:cmax);
    else
        crop = @(M) M;
    end

    %% ---------------------------------------------------------
    %  Von Mises Contour Plot (NO TITLE)
    %% ---------------------------------------------------------
    fig = figure('Name', sprintf('%s 2D Von Mises Stress', datasetName));
    contourf(crop(sigma_vm), 15, 'LineStyle', 'none');
    axis equal tight

    inset = get(gca,'TightInset');
    extra = [0.10 0.10 0.25 0.10];
    set(gca,'LooseInset', inset + extra)

    cb = colorbar;
    cb.Label.String      = 'von Mises stress [MPa]';
    cb.Label.Interpreter = 'none';
    cb.Label.FontSize    = 20;
    colormap(jet);

    set(fig,'PaperPositionMode','auto');
    exportgraphics(fig, ...
        fullfile(saveFolder, sprintf('%s_VM.pdf', datasetName)), ...
        'ContentType', 'vector');

    %% ---------------------------------------------------------
    %  Equivalent 2D Strain
    %% ---------------------------------------------------------
    eps_eq = sqrt( (2/3) * ( (exx - eyy).^2 + exx.^2 + eyy.^2 + 3 * exy.^2 ) );

    fig = figure('Name', sprintf('%s Equivalent 2D Strain', datasetName));
    contourf(crop(eps_eq), 15, 'LineStyle', 'none');
    axis equal tight

    inset = get(gca,'TightInset');
    extra = [0.10 0.10 0.25 0.10];
    set(gca,'LooseInset', inset + extra)

    cb = colorbar;
    cb.Label.String      = "Equivalent strain [-]";
    cb.Label.Interpreter = 'none';
    cb.Label.FontSize    = 25;
    colormap(jet);

    set(fig,'PaperPositionMode','auto');
    exportgraphics(fig, ...
        fullfile(saveFolder, sprintf('%s_eqstrain.pdf', datasetName)), ...
        'ContentType', 'vector');

    %% ---------------------------------------------------------
    %  Plot individual stress and strain fields (NO TITLES)
    %% ---------------------------------------------------------
    for compPlot = 1:length(components)

        comp = components{compPlot};

        switch comp
            case 'exx', strainMatrix = exx;
            case 'eyy', strainMatrix = eyy;
            case 'exy', strainMatrix = exy;
        end

        % -------------------------
        % Stresses for individuals
        % -------------------------
        switch comp
            case 'exx'
                stressMatrix = C * (strainMatrix + nu * eyy);
                stressSymbol = "σxx";
                strainSymbol = "εxx";

            case 'eyy'
                stressMatrix = C * (strainMatrix + nu * exx);
                stressSymbol = "σyy";
                strainSymbol = "εyy";

            case 'exy'
                % *** Corrected shear stress ***
                stressMatrix = 2 * G_MPa .* strainMatrix;
                stressSymbol = "τxy";
                strainSymbol = "εxy";
        end

        %% Stress contour
        fig = figure('Name', sprintf('%s Stress (%s)', datasetName, comp));
        contourf(crop(stressMatrix), 15, 'LineStyle', 'none');
        axis equal tight

        inset = get(gca,'TightInset');
        extra = [0.10 0.10 0.25 0.10];
        set(gca,'LooseInset', inset + extra)

        cb = colorbar;
        cb.Label.String      = stressSymbol + " [MPa]";
        cb.Label.Interpreter = 'none';
        cb.Label.FontSize    = 25;
        colormap(jet);

        set(fig,'PaperPositionMode','auto');
        exportgraphics(fig, ...
            fullfile(saveFolder, sprintf('%s_%s_stress.pdf', datasetName, comp)), ...
            'ContentType', 'vector');

        %% Strain contour
        fig = figure('Name', sprintf('%s Strain (%s)', datasetName, comp));
        contourf(crop(strainMatrix), 15, 'LineStyle', 'none');
        axis equal tight

        inset = get(gca,'TightInset');
        extra = [0.10 0.10 0.25 0.10];
        set(gca,'LooseInset', inset + extra)

        cb = colorbar;
        cb.Label.String      = strainSymbol + " [-]";
        cb.Label.Interpreter = 'none';
        cb.Label.FontSize    = 25;
        colormap(jet);

        set(fig,'PaperPositionMode','auto');
        exportgraphics(fig, ...
            fullfile(saveFolder, sprintf('%s_%s_strain.pdf', datasetName, comp)), ...
            'ContentType', 'vector');
    end
end
