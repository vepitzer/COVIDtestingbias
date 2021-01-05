function sdate = datetime0(idate)
%DATETIME0 Shortcut to datetime when date are given as integer vector

    
    validateattributes(idate, {'numeric'},...
        {'integer','vector','positive','numel',3});
    if idate(2) > 12
        error('Invalid month %g > 12',idate(2))
    end
    if idate(3) > 31
        error('Invalid day %g > 12',idate(3))
    end
    sdate = datetime(datestr(datenum(idate)));

end

