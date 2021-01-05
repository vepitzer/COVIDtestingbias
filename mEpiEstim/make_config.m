function config = make_config(  varargin)
%MAKE_CONFIG Set and check parameter settings of estimate_R

%*************************************************************************%
%   This is MATLAB translation of R function make_config from             %
%   EpiEstim 2.2-3 package                                                % 
%   https://cran.r-project.org/web/packages/EpiEstim/index.html           %
%*************************************************************************%

% 2020/07/25    Created
% 2020/08/04    Add wind_len option

    % Set defaults
    config = [];
    incid = [];
    
    % Intercept object
    if ~isempty(varargin)
        for n = 1:2:length(varargin)
              switch lower(varargin{n})
                  case 'config'
                       config = varargin{n+1};  
                       if ~isa(config,'estimate_R_config')
                           error('Expect estimate_R_config object.')
                       end
                       varargin(n:n+1) = [];
                       break
              end                      
        end
    end
    
    if isempty(config)
        % Create object
        config = estimate_R_config();
    end

    % Scan input
    if ~isempty(varargin)
        for n = 1:2:length(varargin)
            switch lower(varargin{n})
                case 'incid'
                   incid = varargin{n+1};                    
                case 't_start'
                    config.t_start = varargin{n+1};
                    validateattributes(config.t_start, {'numeric'},...
                        {'integer','vector','nonnegative','increasing'});
                case 't_end'
                    config.t_end = varargin{n+1};
                    validateattributes(config.t_end, {'numeric'},...
                        {'integer','vector','nonnegative','increasing'});
                 case 'n1'
                   config.n1 = varargin{n+1};
                    validateattributes(config.n1, {'numeric'},...
                        {'integer','scalar','positive'});
                 case 'n2'
                    config.n2 = varargin{n+1};
                    validateattributes(config.n2, {'numeric'},...
                        {'integer','scalar','positive'});
                case 'mean_si'
                    config.mean_si = varargin{n+1};
                    validateattributes(config.mean_si, {'numeric'},...
                        {'real','scalar','positive'});
                case 'max_mean_si'
                    config.max_mean_si = varargin{n+1};
                    validateattributes(config.max_mean_si, {'numeric'},...
                        {'real','scalar','positive'});  
                case 'max_std_si'
                    config.max_std_si = varargin{n+1};
                    validateattributes(config.max_std_si, {'numeric'},...
                        {'real','scalar','positive'});                     
                case 'min_mean_si'
                    config.min_mean_si = varargin{n+1};
                    validateattributes(config.min_mean_si, {'numeric'},...
                        {'real','scalar','positive'});        
                case 'min_std_si'
                    config.min_std_si = varargin{n+1};
                    validateattributes(config.min_std_si, {'numeric'},...
                        {'real','scalar','positive'});                       
                case 'std_si'
                    config.std_si = varargin{n+1};
                    validateattributes(config.std_si, {'numeric'},...
                        {'real','scalar','positive'});
                case 'std_mean_si'
                    config.std_mean_si = varargin{n+1};
                    validateattributes(config.std_mean_si, {'numeric'},...
                        {'real','scalar','positive'});
                case 'std_std_si'
                    config.std_std_si = varargin{n+1};
                    validateattributes(config.std_std_si, {'numeric'},...
                        {'real','scalar','positive'});                    
                case {'distr_si','si_distr'}
                    config.si_distr = varargin{n+1};
                    validateattributes(config.si_distr, {'numeric'},...
                        {'real','vector','>=',0}); 
                case 'method'
                    % dummy
                    method = varargin{n+1};
                case 'n_sim'
                    % MB
                    config.n_sim = varargin{n+1};
                    validateattributes(config.n_sim, {'numeric'},...
                        {'integer','scalar','nonnegative'});
                case 'wind_len'
                    config.wind_len = varargin{n+1};
                otherwise
                    error('Unknown property %s.',varargin{n})
            end
        end
    end

    % checking and processing incid
    if ~isempty(incid)
        incid = process_I(incid);
        T = length(incid.I);
        

        % filling in / checking t_start and t_end
        if isempty(config.t_start) || isempty(config.t_end)
%             msg = ['Default config will estimate R on weekly sliding windows.',...
%             'To change this change the t_start and t_end arguments.'];
%             fprintf('%s\n',msg)
            config.t_start = 2:T-(config.wind_len - 1);
            config.t_end = config.t_start  +  (config.wind_len - 1);
        else
            check_times(config.t_start, config.t_end, T)
            if ~isscalar(unique(config.t_end - config.t_start))
                wind_len = [];
            end
        end
    end

end

