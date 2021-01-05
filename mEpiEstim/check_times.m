function check_times(t_start, t_end, T)
%CHECK_TIMES Check start end end time vectors
% This only produces errors, does not return anything.
%
%=========================================================================%
% This is MATLAB translation of R function check_times from 
% EpiEstim package https://cran.r-project.org/web/packages/EpiEstim/index.html
%=========================================================================%

% 2020/07/25    Created

    if length(t_start) ~= length(t_end)
        error('t_start and t_end must have the same length.')
    end
    if any(t_start > t_end)
        error('t_start(i) must be <= t_end(i) for all i.')
    end
    if any(t_start < 2) || any(t_start > T )
        error(['t_start must be a vector of integers between 2 ',...
            'and the number of timesteps in incid.'])
    end
    if any(t_end < 2) || any(t_end > T)
        error(['t_end must be a vector of integers between 2',...
            'and the number of timesteps in incid.'])
    end
end
