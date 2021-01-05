function lambda = overall_infectivity(incid, si_distr)
%OVERALL_INFECTIVITY Overall Infectivity Due To Previously Infected 
% Individuals
%
%=========================================================================%
% This is MATLAB translation of R function overall_infectivity from 
% EpiEstim package https://cran.r-project.org/web/packages/EpiEstim/index.html
%=========================================================================%
%
% Input:
%   incid
%       One of the following:
%       - A vector of non-negative integers containing an 
%         incidence time series
%       - A structure of non-negative integers with two columns, so
%         that incid.local} contains the incidence of cases due to local
%         transmission and incid.imported} contains the incidence of imported
%         with incid.local + incid.imported the total incidence) 
%         Note that the cases from the first time step are always all assumed to be
%         imported cases.
%   si_distr
%       Vector of probabilities giving the discrete distribution of
%       the serial interval.
%
% Output:
%   lambda
%       A vector which contains the overall infectivity lambda_t at
%       each time step

% 2020/07/25    Created

    incid = process_I(incid);
    T = length(incid.I);

    check_si_distr(si_distr, "warning")

    m = length(si_distr);
    
    % MB
    if m < T
        si_distr(m+1:T)=0;
    end

    lambda = NaN(T,1);
    for  t = 2:T
        lambda(t) = sum(si_distr(1:t) .* incid.I(t:-1:1),'omitnan');
    end
    
end
