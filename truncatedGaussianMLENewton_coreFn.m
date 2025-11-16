function [mu_est, sigma_est] = truncatedGaussianMLENewton_coreFn(x_data, xmin, tol, max_iter)
    % Truncated Gaussian MLE estimation using Newton's Method
    % Inputs:
    %   x_data - Observed truncated Gaussian data (x >= xmin)
    %   xmin   - Truncation point
    %   tol    - Convergence tolerance (default: 1e-6)
    %   max_iter - Maximum iterations (default: 100)
    % Outputs:
    %   mu_est - Estimated mean
    %   sigma_est - Estimated standard deviation

    if nargin < 3, tol = 1e-11; end
    if nargin < 4, max_iter = 2e5; end

    % Initial guesses
    mu = mean(x_data);
    sigma = std(x_data);
    n = length(x_data);

    for iter = 1:max_iter
        z = (xmin - mu) / sigma;
        phi_z = normpdf(z);
        Phi_z = normcdf(z);

        % Compute function values
        f_mu = mu - (mean(x_data) - sigma * (phi_z / (1 - Phi_z)));
        f_sigma = sigma^2 - (mean((x_data - mu).^2) - sigma * (phi_z / (1 - Phi_z)) * (xmin - mu));

        % Compute Jacobian matrix elements
        dF_mu_mu = 1 + (phi_z / (sigma * (1 - Phi_z)));
        dF_mu_sigma = phi_z / (1 - Phi_z);
        dF_sigma_mu = -2 * sigma * dF_mu_mu;
        dF_sigma_sigma = 2 * sigma - (phi_z / (1 - Phi_z)) * (xmin - mu);

        % Jacobian matrix
        J = [dF_mu_mu, dF_mu_sigma; dF_sigma_mu, dF_sigma_sigma];

        % Newton step: solve J * delta = -F
        delta = -J \ [f_mu; f_sigma];

        % Update values
        mu = mu + delta(1);
        sigma = sigma + delta(2);

        % Convergence check
        if norm(delta) < tol
            break;
        end
    end

    mu_est = mu;
    sigma_est = sigma;
end
