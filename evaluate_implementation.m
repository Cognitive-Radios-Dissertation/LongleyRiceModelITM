% Comprehensive evaluation of Longley-Rice implementation

fprintf('=== Longley-Rice Implementation Evaluation ===\\n\\n');

% Test 1: Basic functionality with different scenarios
fprintf('1. Basic Functionality Tests:\\n');

% Flat earth test
z_flat = 200 * ones(1, 51);
[L_flat, d_flat] = longley_rice_p2p(z_flat, 100, 970, 52, 2.4, struct());
fprintf('   Flat Earth (5km): %.1f dB, Mode: %d (LOS=1, Diff=2, Scat=3)\\n', L_flat, d_flat.mode);

% Hilly terrain
z_hill = 200 + 50 * sin((0:100:5000)/500);
[L_hill, d_hill] = longley_rice_p2p(z_hill, 100, 970, 52, 2.4, struct());
fprintf('   Hill Terrain (5km): %.1f dB, Mode: %d\\n', L_hill, d_hill.mode);

% Different frequencies
[L_2GHz, d_2GHz] = longley_rice_p2p(z_flat, 100, 2000, 52, 2.4, struct());
fprintf('   Flat Earth @ 2GHz: %.1f dB, Mode: %d\\n', L_2GHz, d_2GHz.mode);

[L_100MHz, d_100MHz] = longley_rice_p2p(z_flat, 100, 100, 52, 2.4, struct());
fprintf('   Flat Earth @ 100MHz: %.1f dB, Mode: %d\\n', L_100MHz, d_100MHz.mode);

% Test 2: Validation checks
fprintf('\\n2. Input Validation Tests:\\n');

try
    longley_rice_p2p(z_flat, 100, 19, 52, 2.4, struct()); % Invalid frequency
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

% Get detailed results
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
    [L_conf, d_conf] = longley_rice_p2p(z_flat, 100, 970, 52, 2.4, opts);
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

fprintf('\\n=== Evaluation Complete ===\\n');