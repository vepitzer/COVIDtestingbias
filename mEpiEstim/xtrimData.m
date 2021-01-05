function incid = trimData(incid,ndays)
%TRIMDATA right trim incidence data 
    if isdatetime(ndays)
        tdate = datenum(ndays);
%         tstart = datenum(incid.start_date);
        tend   = datenum(date);
%         if tdate < tstart || tdate > tend
%             error('Invalid date %s.',ndays)
%         end
        ndays = tend - tdate - 1;
    end
    validateattributes(ndays, {'numeric'},...
        {'integer','vector','>',0,'<',length(incid.I)});
    incid.start_date = datetime(datestr(datenum(incid.start_date) + ...
        length(incid.I) - 1 - ndays));
    incid.I=incid.I(end-ndays:end);
end

