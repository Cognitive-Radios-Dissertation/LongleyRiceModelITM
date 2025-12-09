% Test the restructured Longley-Rice implementation

fprintf('=== Testing Restructured Longley-Rice Implementation ===\\n');

% Add src to path
addpath('src');

% Test 1: Basic functionality
fprintf('\\n1. Basic Functionality Tests:\\n');

% Simple flat earth test
z_flat = 200 * ones(1, 51);
[L_flat, details] = longley_rice_p2p(z_flat, 100, 970, 52, 2.4, struct());
fprintf('   Flat Earth (5km): %.1f dB, Mode: %d\\n', L_flat, details.mode);
fprintf('   Details - A_ref: %.1f dB, A_var: %.1f dB\\n', details.A_ref, details.A_var);

% Hilly terrain
z_hill = 200 + 50 * sin((0:100:5000)/500);
[L_hill, details_hill] = longley_rice_p2p(z_hill, 100, 970, 52, 2.4, struct());
fprintf('   Hill Terrain (5km): %.1f dB, Mode: %d\\n', L_hill, details_hill.mode);

% Different frequencies
[L_2GHz, ~] = longley_rice_p2p(z_flat, 100, 2000, 52, 2.4, struct());
fprintf('   Flat Earth @ 2GHz: %.1f dB\\n', L_2GHz);

% Test 2: Validation checks
fprintf('\\n2. Input Validation Tests:\\n');

try
    longley_rice_p2p(z_flat, 100, 19, 52, 2.4, struct()); % Invalid freq
    fprintf('   ERROR: Should have failed frequency validation\\n');
catch ME
    fprintf('   ✓ Frequency validation works: %s\\n', ME.message);
end

try
    longley_rice_p2p(z_flat, 100, 970, 52, 2.4, struct('eps_r', 150)); % Invalid eps_r
    fprintf('   ERROR: Should have failed eps_r validation\\n');
catch ME
    fprintf('   ✓ Epsilon_r validation works: %s\\n', ME.message);
end

try
    short_z = 200 * ones(1, 5); % Too few points
    longley_rice_p2p(short_z, 100, 970, 52, 2.4, struct());
    fprintf('   ERROR: Should have failed profile length validation\\n');
catch ME
    fprintf('   ✓ Profile length validation works: %s\\n', ME.message);
end

% Test 3: Check internal calculations
fprintf('\\n3. Internal Calculation Tests:\\n');
[L, details] = longley_rice_p2p(z_flat, 100, 970, 52, 2.4, struct());
fprintf('   Free Space Loss: %.1f dB\\n', details.L_bf);
fprintf('   Reference Attenuation: %.1f dB\\n', details.A_ref);
fprintf('   Variability: %.1f dB\\n', details.A_var);
fprintf('   Total Loss: %.1f dB\\n', L);
fprintf('   Horizon distances: [%.0f, %.0f] m\\n', details.prop.d_L(1), details.prop.d_L(2));
fprintf('   Effective heights: [%.1f, %.1f] m\\n', details.prop.h_e(1), details.prop.h_e(2));
fprintf('   Terrain roughness: %.1f m\\n', details.prop.delta_h);

% Test 4: Different confidence levels
fprintf('\\n4. Confidence Level Tests:\\n');
for conf = [0.1, 0.5, 0.9]
    opts = struct('conf', conf);
    [L_conf, ~] = longley_rice_p2p(z_flat, 100, 970, 52, 2.4, opts);
    fprintf('   Confidence %.1f: %.1f dB (variation: %.1f dB)\\n', conf, L_conf, L_conf - L_flat);
end

% Test 5: Distance scaling
fprintf('\\n5. Distance Scaling Test:\\n');
for dist_km = [1, 5, 10, 20]
    n_points = round(dist_km * 1000 / 100) + 1;
    z_test = 200 * ones(1, n_points);
    [L_test, d_test] = longley_rice_p2p(z_test, 100, 970, 52, 2.4, struct());
    fprintf('   %2d km: %.1f dB (Mode: %d)\\n', dist_km, L_test, d_test.mode);
end

% Test 6: Test individual functions
fprintf('\\n6. Individual Function Tests:\\n');

% Test hzns
pfl_test.z = z_flat;
pfl_test.xi = 100;
h_g_test = [52, 2.4];
a_eff_test = 6371000 * 1.33; % Typical k-factor
[dL, the] = hzns(pfl_test, h_g_test, a_eff_test);
fprintf('   hzns - Horizon distances: [%.0f, %.0f] m, Angles: [%.4f, %.4f] rad\\n', ...
    dL(1), dL(2), the(1), the(2));

% Test dlthx
delta_h = dlthx(pfl_test);
fprintf('   dlthx - Terrain roughness: %.1f m\\n', delta_h);

% Test avar
prop_test.d_L = dL;
prop_test.the = the;
prop_test.h_e = [52, 2.4];
prop_test.delta_h = delta_h;
prop_test.a_eff = a_eff_test;

prop_params_test.freq_mhz = 970;
prop_params_test.conf = 0.5;
prop_params_test.clim_code = 5;

[A_var, stats] = avar(0.1, prop_test, prop_params_test);
fprintf('   avar - Variability: %.1f dB (Y_total: %.1f dB)\\n', A_var, stats.Y_total);

fprintf('\\n=== Restructured Implementation Test Complete ===\\n');