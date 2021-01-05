function check_config(config, method) 
%CHECK_CONFIG Checks if the method has all the data its need.

%=========================================================================%
% This is MATLAB translation of R function check_config from 
% EpiEstim package https://cran.r-project.org/web/packages/EpiEstim/index.html
%=========================================================================%

% 2020/07/25    Created

    switch lower(method)
        case 'non_parametric_si'
            check_si_distr(config.si_distr,method)
        case 'parametric_si'
            if isempty(config.mean_si)
                error('method parametric_si requires to specify the config.mean_si argument.')
            end
            if isempty(config.std_si)
                error('method parametric_si requires to specify the config.std_si argument.')
            end
        case 'uncertain_si'
        otherwise
            error('Method %s not yet implemented.',method)            
    end
end


