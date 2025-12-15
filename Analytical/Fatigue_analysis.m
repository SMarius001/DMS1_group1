clear; clc;

%% Geometry and material properties
% Specimen 1
    ai = 13.6;          % initial crack size (mm)
    af = 23.9;         % final crack size   (mm)
    F  = 6000;         % applied force (N)
    sigma_ys = 437;
    sigma_uts = 667;
% Specimen 2
    % ai = 13.1;          % initial crack size (mm)
    % af = 24.3;         % final crack size   (mm)
    % F  = 4000;         % applied force (N)
    % sigma_ys = 461;
    % sigma_uts = 675;

b  = 25;            % plate width (mm)
t  = 8;             % thickness (mm)
A  = b * t;       
R  = 0.1;          
alpha = 1;         

%% Insert Paris parameters
% Format:  m, C
 Paris_params = [
     2.71, 1.0e-08
     2.51, 1.67e-12
     3.56, 2.59e-15
     3.08, 1.13e-13
];

n_paris = size(Paris_params,1);

%% Stress calculations
sigma_max = F / A;       
sigma_min = sigma_max * R;
Delta_sigma = sigma_max - sigma_min;
sigma_o = (sigma_ys + sigma_uts) / 2;

%% Newman closure coefficients
A0 = (0.825 - 0.34*alpha + 0.05*alpha^2) * (cos(pi*sigma_max/(2*sigma_o)))^(1/alpha);
A1 = (0.415 - 0.071*alpha) * (sigma_max / sigma_o);
A3 = 2*A0 + A1 - 1;
A2 = 1 - A0 - A1 - A3;

phi = A0 + A1*R + A2*R^2 + A3*R^3;

%% Geometry correction function
f = @(a) sqrt(2*b.*tan(pi*a./(2*b))./(pi*a)) .* ...
         (0.752 + 2.02*(a/b) + 0.37*(1 - sin(pi*a/(2*b))).^3) ./ ...
          cos(pi*a/(2*b));

%% SIF functions
Kmax = @(a) sigma_max .* sqrt(pi.*a) .* f(a);
Kmin = @(a) sigma_min .* sqrt(pi.*a) .* f(a);
DeltaK = @(a) Kmax(a) - Kmin(a);
DeltaKeff = @(a) ((1 - phi) / (1 - R)) .* DeltaK(a);

%% Range for plotting
a_range = linspace(ai, min(b - 1e-3, af), 200);

%% Preallocate results
N_values = zeros(3,1);
N_curves = zeros(3,length(a_range));

%% Loop through the n Paris parameters
for k = 1:n_paris
    m = Paris_params(k,1);
    C = Paris_params(k,2);

    % Paris integrand
    dNda = @(a) 1 ./ (C .* (DeltaKeff(a)).^m);

    % Total cycles
    N_values(k) = integral(dNda, ai, af);

    % N-a curves
    for i = 1:length(a_range)
        if a_range(i) >= ai
            N_curves(k,i) = integral(dNda, ai, a_range(i));
        end
    end

    fprintf("Set %d: m = %.3f, C = %.3e --> N = %.3e cycles\n", ...
        k, m, C, N_values(k));
end

% 1. Plot N vs a
figure; hold on;
colors = lines(n_paris);

for k = 1:n_paris
    plot(N_curves(k,:), a_range)
end

xlabel('Cycles N')
ylabel('Crack length a [mm]')
legend('Seitl et al','Busari et al','Jesus et al','Xin et al')
grid on


%% 2. Plot ΔKeff vs a
figure; hold on; 
K = DeltaKeff(a_range);
a_range = a_range;
plot(a_range, K, 'LineWidth', 2)
xlabel('Crack length a [mm]')
ylabel('\DeltaK_{eff} [MPa mm^{0.5}]')
grid on

%% Compare the size of the plastic zone 
r_p = (1/pi) .* (DeltaKeff(a_range) ./ sigma_ys).^2; % mm
Linear_assum = a_range ./ 8; % mm (12.5% linear assumption)

idx = find(r_p > Linear_assum, 1, 'first');
if ~isempty(idx)
    fprintf('Plastic zone first exceeds the linear assumption at a = %.4f mm\n', a_range(idx));
else
    fprintf('Plastic zone does not exceed the linear assumption in the sampled range.\n');
end
