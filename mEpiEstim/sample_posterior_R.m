function R_sample = sample_posterior_R( R, n, window) 
%SAMPLE_POSTERIOR_R Sample from the posterior R distribution.
%
%=========================================================================%
% This is MATLAB translation of R function sample_posterior_R from 
% EpiEstim package https://cran.r-project.org/web/packages/EpiEstim/index.html
%=========================================================================%
%
% Input:
%   R      -- structure from the estimate_r function function
%   n      -- an integer specifying the number of samples to be taken from the
%               gamma distribution.
%   window -- an integer (or sequence of integers) specifying the window(s) 
%               from which to estimate R. Defaults to the first window. 
%               If multiple windows  are specified, the resulting samples 
%               will be drawn from several distributions.
% 
% Output:
%   R_sample -- n values of R from the posterior R distribution

% 2020/07/25    Created

%   if inherits(R, "estimate_R")
%     error("R must be generated from the estimate_r function.")
%   end
  
  if isempty(n)
      n = 1000;
  end
  if isempty(window)
      window = 1;
  end

  mu    = R.mean_posterior(window);
  sigma = R.std_posterior(window);
  
  % Gamma distribution needs a shape and scale, which are transformations
  % of the mean (mu) and standard deviation (sigma) 
  
  cv    = sigma ./ mu; % coefficient of variation
  shape = 1 ./ (cv .^ 2);
  scale = mu .* cv .^ 2;
  
  R_sample = rgamma(n, shape, scale);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function r = rgamma( n, a, b)
%RGAMMA Wrapper for Matlab function gamrnd
    r = gamrnd(a,b,n,1);
end

