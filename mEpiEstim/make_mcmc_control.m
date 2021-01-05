function mcmc_control = make_mcmc_control( varargin)
%MAKE_MCMC_CONTROL Creates a list of mcmc control parameters to be used in 
% config.mcmc_control, where config is an argument of the 
% estimate_R function. This is used to configure the MCMC chain used to 
% estimate the serial interval within estimate_R (with method 
% "si_from_data").

%**************************************************************************
%   This is MATLAB translation of R function make_mcmc_control from
%   EpiEstim 2.2-3 package
%   https://cran.r-project.org/web/packages/EpiEstim/index.html
%**************************************************************************

% 2020/07/25    Created

    % Create object
    mcmc_control = estimate_R_mcmc_control();

    % Scan input
    if ~isempty(varargin)
        for n = 1:2:length(varargin)
            switch lower(varargin{n})
                case 'burnin'
                    mcmc_control.burnin = varargin{n+1};
                case 'thin'
                    mcmc_control.thin = varargin{n+1};      
                case 'seed'
                    mcmc_control.seed = varargin{n+1};    
                case 'init_pars'
                    mcmc_control.init_pars = varargin{n+1};
                otherwise
                    error('Unknown property %s.',varargin{n})
            end
        end
    end
end