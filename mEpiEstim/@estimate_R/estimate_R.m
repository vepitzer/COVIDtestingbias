classdef estimate_R < handle
    %Estimated Instantaneous Reproduction Number
    %   estimates the reproduction number of an epidemic, given the
    %   incidence time series and the serial interval distribution.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Based on R function estimate_R from EpiEstim pacgage
    % https://cran.r-project.org/web/packages/EpiEstim/index.html
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Usage
    %   obj = estimate_R(incid,method,config [,si_sample[,si_data]])
    %   obj = estimate_R(incid,method,'config',config [,si_sample[,si_data]])
    %   obj = estimate_R('incid',incid,'method',method,'config',config',...
    %           'si_sample',si_sample,...)
    %
    %   incid
    %       a integer vector or the structure with the field
    %           incid.I  (integer vector) or
    %           incid.local (integer vector)
    %           incid.imported (integer vector)
    %       and optional:
    %           incid.start_date (datetime class) or
    %           incid.dates  (datetime class vector)
    %           incid.name
    %           inciod.population
    
    properties
        I                     % integer row vector
        I_local               % integer row vector
        I_imported            % integer row vector
        dates                 % datetime row vector
        start_date            % datetime (MB)
        config                % class 'estimate_R_config'
        method                % char
        si_distr              % double row vector
        t_start               % integer row vector
        t_end                 % integer row vector
        % MB -------------------------------
        name
        population
        IR                    % incidence rate
        AR                    % attack rate
        N = 1e5;              % standard population size
    end
    
    % Results (double row vectors)
    properties (SetAccess = protected)
        mean_posterior
        std_posterior
        quantile_0025_posterior
        quantile_005_posterior
        quantile_025_posterior
        median_posterior
        quantile_075_posterior,
        quantile_095_posterior
        quantile_0975_posterior
        SI_Moments
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods
        function obj = estimate_R(varargin) 
            %ESTIMATE_R Construct an instance of this class
                        
            % minimal check
            narginchk(3,inf)
            
            % default values
            incid = [];
            %si_data = [];
            si_sample = [];
            config = [];
                
            if ~ischar(varargin{1})
                % Fixed order inoput
                incid = process_I(varargin{1});
                if ischar(varargin{2})
                    method = varargin{2};
                    if strcmpi(varargin{3},'config')
                        varargin(3) = [];
                    end
                    if isa(varargin{3},'estimate_R_config')
                        config = varargin{3};
                    else
                        error('Excepted input to be class estimate_R_config.')
                    end
                else
                    error('Excpected input to be char.')
                end
                if length(varargin) > 3
                    si_data = varargin{4};
                end
                if length(varargin) > 4
                    si_sample = varargin{5};
                end
            else        
                % name-value pairs input
                if ~isempty(varargin)
                    for n = 1:2:length(varargin)
                        switch lower(varargin{n})
                            case 'incid'
                                incid = process_I(varargin{n+1});
                            case 'method'
                                method = varargin{n+1};
                                if ~ischar(method)
                                    error('Expect input is char.')
                                end
                            case 'si_data'
                                % si_data = varargin{n+1};
                            case 'si_sample'
                                si_sample = varargin{n+1};
                            case 'config'
                                config = varargin{n+1};
                                if ~isa(config,'estimate_R_config')
                                    error('Expect estimate_R_config object.')
                                end
                            otherwise
                                error('Unknown property %s.',varargin{n})
                        end
                    end
                end
            end
            switch lower(method)
                case {'non_parametric_si', 'parametric_si','uncertain_si'}
                case {'si_from_data','si_from_sample'}
                    error('Method %s not yet implemented.',method)
                otherwise
                    error('Unknown method %s.',method);
            end
            
            config = make_config('incid',incid,'method', method,...
                'config',config);
            %config = process_config(config);
            check_config(config, method)

            %========================================================
            estimate_R_func( obj, incid, method, si_sample , config);
            %========================================================
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = private)
        
        function  estimate_R_func( results, incid, method, si_sample, config)
            %ESTIMATE_R_FUNC Doing the heavy work in estimate_R 
            
            incid = process_I(incid);
            T = length(incid.I);
            
            a_prior = (config.mean_prior / config.std_prior)^2;
            b_prior = config.std_prior^2 / config.mean_prior;
            
            check_times(config.t_start, config.t_end, T)
            nb_time_periods = length(config.t_start);
            
            min_nb_cases_per_time_period = ceil(1 / config.cv_posterior^2 - a_prior);
            incidence_per_time_step  = calc_incidence_per_time_step(...
                incid, config.t_start, config.t_end);
            
            if incidence_per_time_step(1) < min_nb_cases_per_time_period
                msg = ['You''re estimating R too early in the epidemic ',...
                    'to get the desired  posterior CV.'];
                warning('%s Minimal cases per time period is %g',...
                    msg, min_nb_cases_per_time_period);
            end
            
            switch lower(method)
                case 'non_parametric_si'
                    si_uncertainty = 'N';
                    parametric_si = 'N';
                case 'parametric_si'
                    si_uncertainty = 'N';
                    parametric_si = 'Y';
                case 'uncertain_si'
                    si_uncertainty = 'Y';
                    parametric_si = 'Y';
            end
            
            if  strcmpi(si_uncertainty,'Y')
                if strcmpi(parametric_si,'Y')
                    mean_si_sample = -1*ones(1,config.n1);
                    std_si_sample = -1*ones(1,config.n1);
                    for k = 1:config.n1
                        while mean_si_sample(k) < config.min_mean_si || ...
                                mean_si_sample(k) > config.max_mean_si
                            mean_si_sample(k) = rnorm(1, ...
                                config.mean_si,config.std_mean_si);
                        end
                        while std_si_sample(k) < config.min_std_si || ...
                                std_si_sample(k) > config.max_std_si
                            std_si_sample(k) = rnorm(1, ...
                                config.std_si, config.std_std_si);
                        end
                    end
                    temp  = NaN(nb_time_periods,config.n2,config.n1);
                    temp1 = NaN(T,config.n1);
                    for k = 1:config.n1
                        [temp(:,:,k),temp1(:,k)] = sample_from_posterior(config.n2,...
                            incid, mean_si_sample(k), std_si_sample(k),...
                            [], a_prior,b_prior, config.t_start, config.t_end);
                    end
                    config.si_distr = temp1';
                    r_sample = NaN(config.n2 * config.n1, nb_time_periods);
                    for k = 1:config.n1
                        jj = find( config.t_end > mean_si_sample(k));
                        r_sample(((k - 1) * config.n2 + 1):(k * config.n2),jj) = ...
                            temp(jj,:,k)';
                    end
                    results.mean_posterior = mean(r_sample,'omitnan');
                    results.std_posterior  = std( r_sample,'omitnan');
                    results.quantile_0025_posterior = quantile(r_sample, 0.025);
                    results.quantile_005_posterior  = quantile(r_sample, 0.05);
                    results.quantile_025_posterior  = quantile(r_sample, 0.25);
                    results.median_posterior        = quantile(r_sample, 0.5);
                    results.quantile_075_posterior  = quantile(r_sample, 0.75);
                    results.quantile_095_posterior  = quantile(r_sample, 0.95);
                    results.quantile_0975_posterior = quantile(r_sample, 0.975);
                else
                end
            else
                % CertainSI
                if strcmp(parametric_si,'Y')
                    config.si_distr = discr_si(0:T - 1, ...
                        config.mean_si, config.std_si);
                end
                
                if length(config.si_distr) < T + 1
                    config.si_distr(length(config.si_distr) + 1:T + 1) = 0;
                end
                final_mean_si = sum(config.si_distr .* ...
                    (0:length(config.si_distr) -  1));
                Finalstd_si = sqrt(sum(config.si_distr .* ...
                    (0:length(config.si_distr) -  1).^2) - final_mean_si^2);
                post = posterior_from_si_distr( incid, ...
                    config.si_distr, a_prior, b_prior, config.t_start, config.t_end);
                
                a_posterior = post(1,:);
                b_posterior = post(2,:);
                results.mean_posterior = a_posterior .* b_posterior;
                results.std_posterior = sqrt(a_posterior) .* b_posterior;
                
                results.quantile_0025_posterior = ...
                    qgamma(0.025,a_posterior, b_posterior);
                results.quantile_005_posterior  = ...
                    qgamma(0.050,a_posterior, b_posterior);
                results.quantile_025_posterior  = ...
                    qgamma(0.250,a_posterior, b_posterior);
                results.median_posterior  = ...
                    qgamma(0.500,a_posterior, b_posterior);
                results.quantile_075_posterior  = ...
                    qgamma(0.750,a_posterior, b_posterior);
                results.quantile_095_posterior  = ...
                    qgamma(0.950,a_posterior, b_posterior);
                results.quantile_0975_posterior = ...
                    qgamma(0.975,a_posterior, b_posterior);
            end
            results.config = config;
            %             results.mean_posterior = mean_posterior;
            %             results.std_posterior = std_posterior;
            %             results.quantile_0025_posterior = quantile_0025_posterior;
            %             results.quantile_005_posterior = quantile_005_posterior;
            %             results.median_posterior = median_posterior;
            %             results.quantile_095_posterior = quantile_095_posterior;
            %             results.quantile_0975_posterior = quantile_0975_posterior;
            
            results.t_start  = config.t_start;
            results.t_end    = config.t_end;
            results.method   = method;
            results.si_distr = config.si_distr;
            if strcmpi(si_uncertainty,'Y')
                results.SI_Moments = [mean_si_sample;std_si_sample];
            else
                results.SI_Moments = [final_mean_si, Finalstd_si];
            end
%             if isfield(incid,'dates')
%                 results.dates = check_dates(incid);
%                 if isdatetime(results.dates)
%                     % MB
%                     results.start_date = result.dates(1);
%                 end
            if isfield(incid,'start_date')
                % MB
                if ~isempty(incid.start_date)
                    if ~isdatetime(incid.start_date)
                        error('incid.start_date must be an object of class datetime.')
                    end
                    results.start_date = incid.start_date;
                    results.dates =  ...
                        datetime(datestr(datenum(results.start_date): ...
                        (datenum(results.start_date) + T - 1)))';
                end
            end
            if isempty(results.dates)
                results.dates = 1:T;
            end
            results.I = incid.I;
            results.I_local = incid.local;
            results.I_imported = incid.imported;
            
            % additional properties (MB)            
            results.name = incid.name;
            results.population = incid.population;
            if ~isempty(results.population) && ~isempty(config.wind_len)
                nn = length(results.t_start);
                results.IR = NaN(1,nn);
                for n = 1:nn
                    results.IR(n) = sum(results.I(results.t_start(n):results.t_end(n)))/...
                        results.population*results.N;
                end
                results.AR = cumsum(results.I)/results.population*results.N;
                nar = length(results.AR) - length(results.t_start)  + 1;
                results.AR = results.AR(nar:end);
            end
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    methods (Access = public)
        plot(obj,varargin)
        [T,T1] = export(obj,varargin)
        T = exportRes(obj,varargin)
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function incidence_per_time_step = calc_incidence_per_time_step(...
    incid, t_start, t_end)
%CALC_INCIDENCE_PER_TIME_STEP Calculates the cumulative incidence over 
%time steps
    nb_time_periods = length(t_start);
    incidence_per_time_step = zeros(length(nb_time_periods),1);
    for i = 1:length(nb_time_periods)
        incidence_per_time_step(i) = sum(incid.I(t_start(i):t_end(i)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function post = posterior_from_si_distr(incid, si_distr, ...
    a_prior, b_prior, t_start, t_end)
%POSTERIOR_FROM_SI_DISTR Calculates the parameters of the Gamma posterior 
%distribution from the discrete SI distribution

    nb_time_periods = length(t_start);
    lambda = overall_infectivity(incid, si_distr);
    final_mean_si = sum(si_distr .* (0:length(si_distr) - 1));
    a_posterior = NaN(length(nb_time_periods),1);
    b_posterior = NaN(length(nb_time_periods),1);
    for t = 1: nb_time_periods
        if t_end(t) >  final_mean_si
            a_posterior(t) = a_prior + sum(incid.local(t_start(t):t_end(t)));
            % only counting local cases on the 'numerator'
        end
    end
    for t = 1:nb_time_periods
        if t_end(t) >  final_mean_si
            b_posterior(t) = 1 / ...
                (1 / b_prior + sum(lambda(t_start(t):t_end(t))));
        end
    end
    post(1,:) =  a_posterior;
    post(2,:) =  b_posterior;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [sample_r_posterior, si_distr] = sample_from_posterior(...
    sample_size, incid, mean_si, std_si, si_distr , a_prior, b_prior, ...
    t_start, t_end)
%SAMPLE_FROM_POSTERIOR Samples from the Gamma posterior distribution for a
% given mean SI and std SI

    nb_time_periods = length(t_start);
    T = length(incid.I);

    if isempty(si_distr)
        si_distr = discr_si(0:(T - 1), mean_si, std_si);
    end

    final_mean_si =sum(si_distr .* (0: length(si_distr) - 1));
    lambda = overall_infectivity(incid, si_distr);
    a_posterior = NaN(1,nb_time_periods);
    b_posterior = NaN(1,nb_time_periods);
    for t = 1:nb_time_periods
        if t_end(t) > final_mean_si
            a_posterior(t) = a_prior + sum(incid.local(t_start(t):t_end(t)));
            % only counting local cases on the "numerator"
        end
    end
    for t = 1:nb_time_periods
        if t_end(t) >  final_mean_si
            b_posterior(t) = 1 / ...
                (1 / b_prior + sum(lambda(t_start(t):t_end(t)),'omitnan'));
        end
    end
    sample_r_posterior = NaN(nb_time_periods,sample_size);
    for t = 1:nb_time_periods
        if ~isnan(a_posterior(t))
            sample_r_posterior(t,:) = ...
                rgamma(sample_size,a_posterior(t),b_posterior(t));
        end
    end
    if sample_size == 1
        %sample_r_posterior = reshape(sample_r_posterior, 1);
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function res = qgamma(p,a, b)
%QGAMMA Wrapper for Matlab function gaminv
    res = gaminv(p,a,b);
end

function r = rgamma( n, a, b)
%RGAMMA Wrapper for Matlab function gamrnd
    r = gamrnd(a,b,n,1);
end

function r = rnorm( n, mean, std)
%RNORM Wrapper for Matlab function normrnd
    r = normrnd(mean,std,1,n);
end




