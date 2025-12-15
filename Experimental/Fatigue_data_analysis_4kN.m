clc; clear;

%% --- Parameters ---------------------------------------------------------
% Specimen 2
ai = 13.1;         % initial crack length [mm]
af = 24.3;         % final crack size   [mm]
F  = 4000;         % applied force [N]
sigma_ys = 461;    % [MPa]
sigma_uts = 675;   % [MPa]
E = 205000;        % [MPa]
b  = 25;           % plate width [mm]
t  = 8;            % thickness [mm]
A  = b * t;        % cross-section [mm^2]
sigma_max = F / A; % [MPa]
R  = 0.1;          
alpha = 1;   
a0 = 13.1;

%% --- Analytical CMOD function (Tada-based) ------------------------------
V = @(a) (1.46 + 3.42*(1 - cos((pi*a)/(2*b)))) ./ (cos((pi*a)/(2*b))).^2;
delta = @(a) ((4*sigma_max.*a)./E) .* V(a);   % [mm]

% analytical a-range (for plotting)
a_analytic = linspace(13.1, 23.5, 2000);
disp_analytic = delta(a_analytic);

%% --- FEM data -----------------------------------------------------------
displacement_FEM = [2.81E-02
2.91E-02
3.25E-02
3.70E-02
4.23E-02
4.85E-02
5.62E-02
6.55E-02
7.65E-02
8.97E-02
1.07E-01
1.27E-01
1.53E-01
1.86E-01
2.30E-01
2.88E-01
3.66E-01
4.76E-01
6.29E-01
8.66E-01
1.24E+00
1.89E+00
3.02E+00
5.13E+00
];

a_FEM = 13 + [9.2928e-002
0.22147
0.68535
1.1808
1.6617
2.129
2.6211
3.0921
3.5607
4.0476
4.5125
4.9988
5.4942
5.9684
6.4604
6.9227
7.3847
7.8719
8.333
8.821
9.2969
9.7652
10.196
10.616
];

DeltaK_FEM = [103.4
181.84
220.27
231.41
249.67
264.28
288.74
311.85
334.62
359.91
383.34
412.67
444.11
478.18
510.86
548.87
593.98
625.34
664.07
689.34
739.8
779.88
862.41
935.62
];

% Verify FEM arrays are same length
if numel(displacement_FEM) ~= numel(a_FEM) || numel(a_FEM) ~= numel(DeltaK_FEM)
    error('FEM arrays displacement_FEM, a_FEM and DeltaK_FEM must be same length');
end

% Ensure monotonicity for interpolation; if not, sort them
if ~all(diff(displacement_FEM) > 0)
    [displacement_FEM, idx] = sort(displacement_FEM);
    a_FEM = a_FEM(idx);
    DeltaK_FEM = DeltaK_FEM(idx);
    fprintf('Note: displacement_FEM was not monotonic; data sorted for interpolation\n');
end

%% --- Experimental CMOD --------------------------------------------------

displacement_exp = [
0.069138672
0.068606738
0.067915808
0.06842402
0.068298082
0.068181757
0.069114145
0.06911483
0.068578448
0.069625925
0.068993728
0.06885096
0.068304803
0.068921786
0.068892564
0.068758692
0.068596583
0.068997752
0.069231312
0.069101285
0.068901565
0.069332491
0.069419075
0.069444339
0.06892005
0.0683357
0.069134038
0.069154184
0.068656053
0.068353075
0.068574864
0.069507513
0.068637434
0.069153655
0.069299903
0.069282014
0.069240365
0.068257447
0.068203922
0.069093983
0.068587679
0.068494264
0.069904919
0.069528889
0.069593157
0.068877514
0.069388535
0.069747042
0.070502024
0.06906617
0.06924406
0.069201681
0.069214042
0.069937419
0.070040282
0.07014272
0.069538165
0.069391411
0.070577469
0.06921931
0.06990521
0.069638539
0.069926362
0.070027139
0.070083897
0.07067227
0.0699188
0.069915678
0.070052836
0.069794576
0.070041459
0.070863414
0.070273775
0.071058493
0.071081582
0.070459489
0.070227939
0.071135674
0.071420852
0.070603769
0.071222219
0.070822026
0.071520392
0.071179803
0.070657465
0.070383657
0.072864983
0.071819294
0.071232799
0.071650598
0.072125669
0.071807731
0.071967687
0.072099667
0.07261559
0.07256763
0.071673307
0.072188873
0.07295673
0.073793981
0.073402923
0.072945874
0.072755258
0.073213179
0.073457424
0.074211474
0.073970068
0.07455897
0.076250907
0.075799871
0.075175066
0.076602314
0.076373372
0.076332707
0.076880995
0.077559758
0.077141766
0.078187783
0.078714956
0.079866488
0.08050875
0.08128823
0.081379455
0.082866732
0.083833177
0.084803898
0.086079222
0.087425892
0.089050923
0.090958779
0.093195555
0.096575988
0.097947267
0.100883832
0.104653051
0.108865669
0.114811985
0.120473681
0.130066982
0.137906126
0.150998749
0.169655421
0.192868272
0.234668341
0.299948886
0.459749573
0.802845459
1.259329144
1.727852789
2.143767941
2.449417475
3.429861717];

N_exp = (1000:1000:152000)'; % vector matches length
if numel(N_exp) ~= numel(displacement_exp)
    error('Length mismatch between N_exp and displacement_exp');
end

%% --- Smooth noisy CMOD for stability ------------------------------------
displacement_exp_s = smoothdata(displacement_exp,'loess',7);
nExp = numel(displacement_exp_s);

% --- Allocate crack length estimates -------------------------------------
a_from_analytic = nan(nExp,1);
a_from_fem      = nan(nExp,1);

%% --- Calibration of experimental CMOD -----------------------------------
% Baseline offset using analytic value at a0
d_ref_a0 = delta(a0);         % analytic delta at a0
% Mean of first few stable experimental points - keep original choice but warn
nBaselinePts = 20;
if nBaselinePts >= nExp
    error('nBaselinePts is too large for displacement_exp length');
end
baseline_exp = mean(displacement_exp_s(1:nBaselinePts));
offset = baseline_exp - d_ref_a0;

% correct experimental
d_exp_cal = displacement_exp_s - offset;

%% --- 1) Invert analytic CMOD --------------------------------------------
a_min = 13; 
a_max = 23.5;
opts = optimset('Display','off');

for i = 1:nExp
    d = d_exp_cal(i);
    % If d outside analytic delta range, skip early
    d_min = min(delta(a_analytic));
    d_max = max(delta(a_analytic));
    if d < d_min || d > d_max
        a_from_analytic(i) = NaN;
        continue;
    end
    % Check bracket
    fa = delta(a_min) - d;
    fb = delta(a_max) - d;
    if fa * fb < 0
        % bracketed: use fzero with bracket
        try
            a_sol = fzero(@(a) delta(a) - d, [a_min a_max], opts);
            a_from_analytic(i) = a_sol;
        catch
            a_from_analytic(i) = NaN;
        end
    else
        % Not bracketed: fallback to bounded minimization of abs error
        funmin = @(a) abs(delta(a) - d);
        try
            a_sol = fminbnd(funmin, a_min, a_max, optimset('TolX',1e-6,'Display','off'));
            if abs(delta(a_sol) - d) < 1e-6
                a_from_analytic(i) = a_sol;
            else
                a_from_analytic(i) = NaN;
            end
        catch
            a_from_analytic(i) = NaN;
        end
    end
end

% Report how many mapped
fprintf('Analytic inversion: %d / %d points mapped (non-NaN)\n', sum(isfinite(a_from_analytic)), nExp);

%% --- 2) Invert FEM curve by interpolation -------------------------------
minF = min(displacement_FEM);
maxF = max(displacement_FEM);

for i = 1:nExp
    d = d_exp_cal(i);
    if d >= minF && d <= maxF
        a_from_fem(i) = interp1(displacement_FEM, a_FEM, d, 'pchip');
    else
        a_from_fem(i) = NaN;
    end
end

fprintf('FEM interpolation: %d / %d points mapped (non-NaN)\n', sum(isfinite(a_from_fem)), nExp);

%% --- SIF from experiment ------------------------------------------------
sigma_min = sigma_max * R;
Delta_sigma = sigma_max - sigma_min;
sigma_o = (sigma_ys + sigma_uts) / 2;

% Newman closure coefficients
A0 = (0.825 - 0.34*alpha + 0.05*alpha^2) * (cos(pi*sigma_max/(2*sigma_o)))^(1/alpha);
A1 = (0.415 - 0.071*alpha) * (sigma_max / sigma_o);
A3 = 2*A0 + A1 - 1;
A2 = 1 - A0 - A1 - A3;

phi = A0 + A1*R + A2*R^2 + A3*R^3;

% Geometry correction function
f = @(a) sqrt(2*b.*tan(pi*a./(2*b))./(pi*a)) .* ...
         (0.752 + 2.02*(a/b) + 0.37*(1 - sin(pi*a/(2*b))).^3) ./ ...
          cos(pi*a./(2*b));

% SIF functions
Kmax = @(a) sigma_max .* sqrt(pi.*a) .* f(a);
Kmin = @(a) sigma_min .* sqrt(pi.*a) .* f(a);
DeltaK = @(a) Kmax(a) - Kmin(a);
DeltaKeff = @(a) ((1 - phi) / (1 - R)) .* DeltaK(a);

% compute DeltaK for analytic mapped points
DeltaK_analytic = DeltaKeff(a_from_analytic);

%% --- Compute Paris data for analytic inversion --------------------------
[daN_A, DKmid_A, nDiscardA] = computeParisData(a_from_analytic, N_exp, DeltaK_analytic);
fprintf('Analytic Paris data: %d pairs used, %d discarded due to NaN or invalid.\n', numel(daN_A), nDiscardA);

%% --- Compute Paris data for FEM inversion -------------------------------
DeltaK_FEM_interp = interp1(a_FEM, DeltaK_FEM, a_from_fem, 'pchip', NaN);
[daN_F, DKmid_F, nDiscardF] = computeParisData(a_from_fem, N_exp, DeltaK_FEM_interp);
fprintf('FEM Paris data: %d pairs used, %d discarded due to NaN or invalid.\n', numel(daN_F), nDiscardF);

%% --- Trim endpoints to avoid large DeltaK values ------------------------
nTrim_beg = 50;
nTrim_end = 5;
% Trim emp inversion
daN_A  = daN_A((1+nTrim_beg) : end-nTrim_end);
DKmid_A = DKmid_A((1+nTrim_beg) : end-nTrim_end);

% Trim FEM inversion
daN_F  = daN_F((1+nTrim_beg) : end-nTrim_end);
DKmid_F = DKmid_F((1+nTrim_beg) : end-nTrim_end);

%% --- Paris fits with bootstrap for CI -----------------------------------
nBoot = 1000; % bootstrap draws

[m_A, C_A, m_boot_A, C_boot_A] = parisFitBootstrap(DKmid_A, daN_A, nBoot);
[m_F, C_F, m_boot_F, C_boot_F] = parisFitBootstrap(DKmid_F, daN_F, nBoot);

% parameter CIs (for reporting)
mCI_A = prctile(m_boot_A, [2.5 97.5]);
CCI_A = prctile(C_boot_A, [2.5 97.5]);
mCI_F = prctile(m_boot_F, [2.5 97.5]);
CCI_F = prctile(C_boot_F, [2.5 97.5]);

fprintf('Analytic fit: m = %.3f (%.3f, %.3f), C = %.3e (%.3e, %.3e)\n', ...
    m_A, mCI_A(1), mCI_A(2), C_A, CCI_A(1), CCI_A(2));

fprintf('FEM fit: m = %.3f (%.3f, %.3f), C = %.3e (%.3e, %.3e)\n', ...
    m_F, mCI_F(1), mCI_F(2), C_F, CCI_F(1), CCI_F(2));

%% --- Plot Paris data with bootstrap CI envelope -------------------------
plotParisBootstrap(DKmid_A, daN_A, m_A, C_A, m_boot_A, C_boot_A, ...
   'Paris Fit using \DeltaK\_analytic');

plotParisBootstrap(DKmid_F, daN_F, m_F, C_F, m_boot_F, C_boot_F, ...
   'Paris Fit using \DeltaK\_FEM');

%% --- Plot results -------------------------------------------------------
figure; hold on;
plot(a_analytic, disp_analytic, 'LineWidth', 1.5);                               % analytic curve
maskA = isfinite(a_from_analytic) & isfinite(d_exp_cal);
plot(a_from_analytic(maskA), d_exp_cal(maskA), '.', 'MarkerEdgeColor','blue')    % mapped analytical
plot(a_FEM, displacement_FEM, 'o-');                                             % FEM points
maskF = isfinite(a_from_fem) & isfinite(d_exp_cal);
plot(a_from_fem(maskF), d_exp_cal(maskF), '.', 'MarkerEdgeColor','red');         % mapped FEM
xlabel('Crack length a [mm]');
ylabel('CMOD [mm]');
legend('Analytic \delta(a)', 'Experimental mapped from analytical', 'FEM', 'Experimental mapped from FEM', 'Location','best');
grid on;

% 1. Plot a vs N for all (only finite values)
figure; hold on;
if any(isfinite(a_from_analytic))
    plot(N_exp(isfinite(a_from_analytic)), a_from_analytic(isfinite(a_from_analytic)));
end
if any(isfinite(a_from_fem))
    plot(N_exp(isfinite(a_from_fem)), a_from_fem(isfinite(a_from_fem)));
end
xlabel('Cycles N');
ylabel('Crack length a [mm]');
legend('Mapped to analytical','Mapped to FEM');
grid on;

%% --- Functions ----------------------------------------------------------

function [da_dN, DK_mid, nDiscard] = computeParisData(a_vec, N_vec, DK_vec)
    % compute da/dN and mid-point DK
    % remove entries where a_vec or DK_vec are not finite
    valid = isfinite(a_vec) & isfinite(DK_vec);
    a_vec = a_vec(valid);
    N_vec = N_vec(valid);
    DK_vec = DK_vec(valid);

    % ensure ordering by N; if not monotonic, sort by N
    if ~issorted(N_vec)
        [N_vec, idx] = sort(N_vec);
        a_vec = a_vec(idx);
        DK_vec = DK_vec(idx);
    end

    % compute diffs
    % require at least two points
    if numel(a_vec) < 2
        da_dN = [];
        DK_mid = [];
        nDiscard = 0;
        return;
    end

    da = diff(a_vec);
    dN = diff(N_vec);
    da_dN_all = da ./ dN;

    % DK midpoint corresponding to each da/dN
    DK_mid_all = 0.5 * (DK_vec(1:end-1) + DK_vec(2:end));

    % filter out nonpositive values (report count)
    mask = isfinite(da_dN_all) & isfinite(DK_mid_all) & (da_dN_all > 0) & (DK_mid_all > 0);
    da_dN = da_dN_all(mask);
    DK_mid = DK_mid_all(mask);

    nDiscard = sum(~mask) + sum(~valid);
end

function [m, C, m_boot, C_boot] = parisFitBootstrap(DK, da_dN, nBoot)
    % Fit log10(da/dN) = m*log10(DK) + log10(C)
    % Return fitted m and C, and bootstrap distributions for m and C
    mask = isfinite(DK) & isfinite(da_dN) & DK > 0 & da_dN > 0;
    DK0 = DK(mask);
    da0 = da_dN(mask);

    if numel(DK0) < 2
        m = NaN; C = NaN; m_boot = []; C_boot = [];
        warning('Not enough data points to fit Paris law.');
        return;
    end

    x = log10(DK0(:));
    y = log10(da0(:));

    p = polyfit(x, y, 1);
    m = p(1);
    C = 10^(p(2));

    % bootstrap parameter distributions
    n = numel(x);
    m_boot = zeros(nBoot,1);
    C_boot = zeros(nBoot,1);
    for b = 1:nBoot
        idx = randi(n, n, 1);
        pb = polyfit(x(idx), y(idx), 1);
        m_boot(b) = pb(1);
        C_boot(b) = 10^(pb(2));
    end
end

function plotParisBootstrap(DK, daN, m, C, m_boot, C_boot, ttl)
    % Plot data, fit and bootstrap CI envelope computed from bootstrapped parameters
    mask = DK > 0 & daN > 0 & isfinite(DK) & isfinite(daN);
    DK = DK(mask);
    daN = daN(mask);

    if isempty(DK)
        warning('No valid DK/daN points to plot for %s', ttl);
        return;
    end

    figure; hold on;
    set(gca, 'XScale', 'log', 'YScale', 'log');

    DKf = logspace(log10(min(DK)), log10(max(DK)), 500);

    % Compute CI envelope: for each bootstrap pair predict da/dN at DKf, then percentile
    if ~isempty(m_boot) && ~isempty(C_boot)
        Yboot = zeros(numel(m_boot), numel(DKf));
        for bi = 1:numel(m_boot)
            Yboot(bi, :) = C_boot(bi) .* (DKf .^ m_boot(bi));
        end
        lo = prctile(Yboot, 2.5, 1);
        hi = prctile(Yboot, 97.5, 1);
        fill([DKf fliplr(DKf)], [lo fliplr(hi)], [0.8 0.85 1], 'EdgeColor','none', 'FaceAlpha', 0.5);
    end

    loglog(DK, daN, 'o', 'MarkerFaceColor','b', 'MarkerEdgeColor','k');

    dafit = C * DKf.^m;
    loglog(DKf, dafit, 'k-', 'LineWidth', 2);

    xlabel('\DeltaK [MPa \cdot mm^{0.5}]');
    ylabel('da/dN [mm/cycle]');
    title(ttl);
    grid on;
    legend('CI band','data','fit','Location','best');

    % annotate fit results
    mCI = prctile(m_boot, [2.5 97.5]);
    CCI = prctile(C_boot, [2.5 97.5]);
    txt = {
        sprintf('m = %.3f', m)
        sprintf('C = %.3e', C)
        sprintf('m CI = [%.3f, %.3f]', mCI(1), mCI(2))
        sprintf('C CI = [%.3e, %.3e]', CCI(1), CCI(2))
    };

    annotation('textbox',[0.15 0.65 0.2 0.2], ...
               'String',txt, ...
               'FitBoxToText','on', ...
               'BackgroundColor',[1 1 1 0.7], ...
               'EdgeColor','k', ...
               'FontSize',10);
end