function dates = check_dates(incid)
%CHEK_DATES

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on R function check_dates from EpiEstim pacgage
% https://cran.r-project.org/web/packages/EpiEstim/index.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2020/07/29    Created

    dates = incid.dates;
    if isempty(dates)
        return
    end
    if isnumeric(dates)
        validateattributes(dates, {'numeric'},...
            {'integer','vector','nonnegative'});
        fail = unique(diff(dates)) ~= 1;
    elseif isdatetime(dates)
        fail = unique(diff(datenum(dates))) ~= 1;
    else
        error('incid.dates must be an object of class datetime or numeric.')
    end
    if fail
        error('incid.dates must contain dates which are all in a row.')
    end
    
end
