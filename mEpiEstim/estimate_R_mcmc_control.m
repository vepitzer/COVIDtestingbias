classdef estimate_R_mcmc_control < handle
    %ESTIMATE_R_MCMC_CONTROL MCMC control parameters to be used in
    %config.mcmc_control.

    %=========================================================================%
    % Based on R function mcmc_control from 
    % EpiEstim package https://cran.r-project.org/web/packages/EpiEstim/index.html
    %=========================================================================%

    % burnin
    %   A positive integer giving the burnin used in the MCMC when
    %   estimating the serial interval distribution.
    % thin
    %   A positive integer corresponding to thinning parameter; the MCMC
    %   will be run for burnin+n1*thin iterations} 1 in thin
    %   iterations will be recorded, after the burnin phase, so the posterior
    %   sample size is n1.
    % seed
    %   An integer used as the seed for the random number generator at
    %   the start of the MCMC estimation; useful to get reproducible results.
    % init_pars
    %   vector of size 2 corresponding to the initial values of
    %   parameters to use for the SI distribution. This is the shape and scale for
    %   all but the lognormal distribution, for which it is the meanlog and
    %   sdlog.
    
    % 2020/07/29    Created
    
    properties
        burnin = 3000
        thin = 10
        seed = fix(now) %= as.integer(Sys.time()), 
        init_pars 
    end
    
    methods
        function obj = estimate_R_mcmc_control()
            %ESTIMATE_R_MCMC_CONTROL Construct empty object
        end
    end
end

