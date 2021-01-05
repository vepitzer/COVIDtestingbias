function check_si_distr(si_distr, sumToOne, method) 
%CHECK_SI_DISTR  Checks si_distr input                    
% This only produces warnings and errors, does not return anything

%=========================================================================%
% This is MATLAB translation of R function check_si_distr from 
% EpiEstim package https://cran.r-project.org/web/packages/EpiEstim/index.html
%=========================================================================%

% 2020/07/25    Created

    if nargin < 3
        method = 'non_parametric_si';
    end

    if isempty(si_distr)
        error('si_distr argument is missing but is required for method %s.', ....
            method)
    end
    
    validateattributes(si_distr, {'numeric'},...
        {'real','vector','nonnegative'});

    if si_distr(1) ~= 0
        error('si_distr should be so that si_distr(1) = 0.')
    end

    if abs(sum(si_distr) - 1) > 0.01
        switch lower(sumToOne)
            case 'error'
                error('si_distr must sum to 1.')
            case 'warning'
                warning('si_distr does not sum to 1.')
            otherwise
                error('Invalid sumToOne string %s.',sumToOne)
        end
    end

end

