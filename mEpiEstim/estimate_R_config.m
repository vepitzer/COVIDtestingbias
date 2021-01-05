classdef estimate_R_config < handle
    %ESTIMATE_R_CONFIG Creates an emmpty estimate_R_config object

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Based on R function make_config from EpiEstim pacgage
    % https://cran.r-project.org/web/packages/EpiEstim/index.html
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % t_start -- Vector of positive integers giving the starting times of each
    % window over which the reproduction number will be estimated. These must be in
    %  ascending order, and so that for all {i, t_start(i)<=t_end(i)}.
    %  t_start[1] should be strictly after the first day with non null incidence.
    %
    % t_end --Vector of positive integers giving the ending times of each
    % window over which the reproduction number will be estimated. These must be
    % in ascending order, and so that for all i, t_start(i)<=t_end(i).
    %
    % n1 --For method "uncertain_si" and "si_from_data"; positive integer
    % giving the size of the sample of SI distributions to be drawn (see details).
    %
    % n2 -- For methods "uncertain_si", "si_from_data" and "si_from_sample";
    % positive integer giving the size of the sample drawn from the posterior
    % distribution of R for each serial interval distribution considered (see
    % details).}
    %
    % mean_si -- For method "parametric_si" and "uncertain_si" ; positive real
    % giving the mean serial interval (method "parametric_si") or the average mean
    % serial interval (method "uncertain_si", see details).
    %
    % std_si --For method "parametric_si" and "uncertain_si" ; non negative
    % real giving the standard deviation of the serial interval
    % (method "parametric_si") or the average standard deviation of the serial
    % interval (method "uncertain_si", see details).
    %
    % std_mean_si --{For method "uncertain_si" ; standard deviation of the
    % distribution from which mean serial intervals are drawn (see details).}
    %
    % min_mean_si --{For method "uncertain_si" ; lower bound of the
    % distribution from which mean serial intervals are drawn (see details).}
    %
    % max_mean_si -- {For method "uncertain_si" ; upper bound of the
    % distribution from which mean serial intervals are drawn (see details).}
    %
    % std_std_si -- {For method "uncertain_si" ; standard deviation of the
    % distribution from which standard deviations of the serial interval are drawn
    % (see details).}
    %
    % min_std_si -- {For method "uncertain_si" ; lower bound of the distribution
    %  from which standard deviations of the serial interval are drawn (see
    %  details).}
    %
    %  max_std_si --{For method "uncertain_si" ; upper bound of the distribution
    %  from which standard deviations of the serial interval are drawn (see
    %  details).}
    %
    % si_distr -- For method "non_parametric_si" ; vector of probabilities
    % giving the discrete distribution of the serial interval, starting with
    % si_distr(1) (probability that the serial interval is zero), which
    % should be zero.}
    %
    % si_parametric_distr --{For method "si_from_data" ; the parametric
    % distribution to use when estimating the serial interval from data on dates of
    %  symptoms of pairs of infector/infected individuals (see details).
    % Should be one of "G" (Gamma), "W" (Weibull), "L" (Lognormal), "off1G" (Gamma
    % shifted by 1), "off1W" (Weibull shifted by 1), or "off1L" (Lognormal shifted
    % by 1).}
    %
    % mcmc_control -- {An object of class \code{estimate_R_mcmc_control}, as
    % returned by function \code{make_mcmc_control}. }
    %
    % seed -- {An optional integer used as the seed for the random number
    % generator at the start of the function (then potentially reset within the
    % MCMC for method \code{si_from_data}); useful to get reproducible results.}
    %
    % mean_prior -- {A positive number giving the mean of the common prior
    % distribution for all reproduction numbers (see details).}
    %
    % std_prior --{A positive number giving the standard deviation of the
    % common prior distribution for all reproduction numbers (see details).}
    %
    % cv_posterior --{A positive number giving the aimed posterior coefficient
    % of variati
    
    % NEW
    % wind_len -- length of time windov (7 by default)
    
    % 2020/07/29    Created
    % 2020/08/04    Add wind_length option
    
    properties
        t_start
        t_end
        n1 = 500
        n2 = 50
        mean_si
        std_si
        std_mean_si
        min_mean_si
        max_mean_si
        std_std_si
        min_std_si
        max_std_si
        si_distr
        si_parametric_distr
        mcmc_control
        seed
        mean_prior = 5
        std_prior = 5
        cv_posterior = 0.3
        wind_len = 7;  % sliding window length (7 days by default)
    end
    
    methods
        function obj = estimate_R_config()
            %CONFIG Construct a default object
            obj.mcmc_control = estimate_R_mcmc_control();
        end
    end
    
    methods
        function set.t_start(obj,t_start)
            validateattributes(t_start, {'numeric'},...
                {'integer','vector','nonnegative'});
            obj.t_start = t_start;
        end
        function set.t_end(obj,t_end)
            validateattributes(t_end, {'numeric'},...
                {'integer','vector','nonnegative'});
            obj.t_end = t_end;
        end        
        function set.mean_si(obj,mean_si) 
            validateattributes(mean_si,  {'numeric'},...
                {'real','scalar','>',1,});
            obj.mean_si = mean_si;
        end
        function set.std_si(obj,mean_si) 
            validateattributes(mean_si,  {'numeric'},...
                {'real','scalar','>',0,});
            obj.std_si = mean_si;
        end   
        function set.cv_posterior(obj,cv_posterior) 
            validateattributes(cv_posterior,  {'numeric'},...
                {'real','scalar','>',0,});
            obj.cv_posterior = cv_posterior;
        end    
        function set.mean_prior(obj,mean_prior) 
            validateattributes(mean_prior,  {'numeric'},...
                {'real','scalar','>',0,});
            obj.mean_prior = mean_prior;
        end  
        function set.std_prior(obj,std_prior) 
            validateattributes(std_prior,  {'numeric'},...
                {'real','scalar','>',0,});
            obj.std_prior = std_prior;
        end  
        function set.wind_len(obj,wind_len)
            validateattributes(wind_len,  {'numeric'},...
                {'integer','scalar','>',1,});
            obj.wind_len = wind_len;
        end
    end

end
