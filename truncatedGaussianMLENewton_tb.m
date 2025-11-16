%% Example Usage
% Generate synthetic truncated Gaussian data
rng(42);
x_data_full = normrnd(5, 2, [1e6, 1]); % Full Gaussian sample
xmin = 5; % Truncation point
x_data = x_data_full(x_data_full >= xmin); % Apply truncation

% Solve using Newton's method
[mu_est, sigma_est] = truncatedGaussianMLENewton_coreFn(x_data, xmin);

fprintf('Estimated mu: %.4f\n', mu_est);
fprintf('Estimated sigma: %.4f\n', sigma_est);
