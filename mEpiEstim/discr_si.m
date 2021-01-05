function res = discr_si( k, mu, sigma) 
%DISCR_SI Discretized Generation Time Distribution Assuming A Shifted Gamma 
% Distribution
%
%=========================================================================%
% This is MATLAB translation of R function discr_si from EpiEstim package
% https://cran.r-project.org/web/packages/EpiEstim/index.html
%=========================================================================%
%
% Input:
%   k ----- Positive integer, or vector of positive integers for which the 
%           discrete distribution is desired.
%   mu ---- A positive real giving the mean of the Gamma distribution.
%   sigma - A non-negative real giving the standard deviation of the Gamma 
%           distribution.
%
% Output:
%   res --- Gives the discrete probability w(k} that the serial interval is 
%           equal to {k}.

% 2020/07/25    Created

    if sigma < 0
        error("sigma must be >=0.")
    end
    if mu <= 1
        error("mu must be >1.")
    end
    if any(k < 0)
        error("all values in k must be >=0.")
    end

    a = ((mu - 1) / sigma)^2;  % shape
    b = sigma^2 / (mu - 1);    % scale

    res = k .* cdf_gamma(k, a, b) + ...
        (k - 2) .* cdf_gamma(k - 2, a, b) - 2 * (k - 1) .* cdf_gamma(k - 1, a, b);
    res = res + a * b * (2 * cdf_gamma(k - 1, a + 1, b) - ...
        cdf_gamma(k - 2, a + 1, b) - cdf_gamma(k, a + 1, b));
    res(res < 0) = 0;

end
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = cdf_gamma( k, a, b)
%CDF_GAMMA Wraper for Matlab function gamcdf
    p = gamcdf( k, a, b);
end
