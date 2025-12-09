% Create Mock X.04 file (Distance, Height)
d = 0:100:50000; % 50 km path, 100m step
% Create a "double hill" terrain
z = 200 + 100 * sin(d/5000) + 50 * cos(d/1000) + 20 * randn(size(d));
z(1) = 200; % Fix start
z(end) = 200; % Fix end

data = [d', z'];
writematrix(data, 'X.04', 'Delimiter', ' ');
fprintf('Created mock X.04 file with 50km terrain.\n');
