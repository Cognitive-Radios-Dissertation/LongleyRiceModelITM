% RUN_SIMULATION - Generate Pathloss vs Distance for Longley-Rice

% Add src to path
addpath('../src');

% 1. Load Data
data_dir = '../data';
if ~exist(data_dir, 'dir'), mkdir(data_dir); end
filename = fullfile(data_dir, 'X.04');

if ~isfile(filename)
    disp('X.04 not found in data/. Creating mock data for testing...');
    d = 0:100:50000; % 50 km, 100m step
    z = 200 + 150 * sin(d/8000).^2 + 20 * randn(size(d)); % Some hills
    data = [d', z'];
    writematrix(data, filename, 'FileType', 'text', 'Delimiter', ' ');
else
    disp(['Loading ', filename, '...']);
    data = readmatrix(filename, 'FileType', 'text');
end

% Extract and Validate
d_vec = data(:, 1);
z_vec = data(:, 2);

% Uniformity Check
steps = diff(d_vec);
if std(steps) > 1e-3
    disp('Warning: Input profile has irregular spacing. Resampling to constant step...');
    mean_step = mean(steps);
    d_new = (d_vec(1):mean_step:d_vec(end))';
    z_new = interp1(d_vec, z_vec, d_new, 'linear');
    
    d_vec = d_new;
    z_vec = z_new;
    xi = mean_step;
else
    xi = steps(1);
    disp(['Profile step size confirmed: ', num2str(xi), ' meters.']);
end

% 2. Setup Parameters (User Specified)
freq = 970;  % MHz
h_tx = 52;   % meters
h_rx = 2.4;  % meters
conf = 0.5;  % 50% confidence (Median)

options.pol = 0; % Horizontal
options.conf = conf;
options.clim = 5; % Continental Temperate

% 3. Simulation Loop
targets = d_vec; 

loss_results = nan(size(targets)); 
mode_results = nan(size(targets));

disp(['Running Simulation (Total Points: ', num2str(length(targets)), ')...']);

% Pre-calculate FSPL for all points
c = 299.792458;
lambda = c / freq;
fspl_results = 20 * log10(4 * pi * max(targets, 0.1) / lambda);

for i = 1:length(targets)
    d_target = targets(i);
    
    % FIX: Allow simulation for short paths (i < 10)
    % We process the profile up to index i.
    sub_z = z_vec(1:i);
    
    % Run Model
    try
        [L_b, details] = longley_rice_p2p(sub_z, xi, freq, h_tx, h_rx, options);
        loss_results(i) = L_b;
        if isfield(details, 'mode') && ~isempty(details.mode)
             mode_results(i) = details.mode(end);
        end
    catch ME
        % disp(['Error at d=', num2str(d_target), ': ', ME.message]);
    end
end

% 4. Plot Results
fig = figure('Visible', 'off'); 

% Plot ITM
plot(targets/1000, loss_results, 'b-', 'LineWidth', 2, 'DisplayName', 'Longley-Rice (ITM)');
hold on;

% Plot FSPL
plot(targets/1000, fspl_results, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Free Space Path Loss');

% Formatting
grid on;
legend('Location', 'best');
xlabel('Distance (km)');
ylabel('Path Loss (dB)');
title(['Longley-Rice Path Loss vs Distance', char(10), ...
       'Freq: ', num2str(freq), ' MHz, Tx: ', num2str(h_tx), 'm, Rx: ', num2str(h_rx), 'm']);

% Save Plot
results_dir = '../results';
if ~exist(results_dir, 'dir'), mkdir(results_dir); end
save_path = fullfile(results_dir, 'pathloss_plot.png');
saveas(fig, save_path);
disp(['Simulation Complete. Plot saved to ', save_path]);
