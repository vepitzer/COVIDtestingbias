function res = process_I(incid) 
%PROCESS_I Convert the input 'incid' to a structure that EpiEstim
%understands.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on R function process_I from EpiEstim pacgage
% https://cran.r-project.org/web/packages/EpiEstim/index.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Input:
%   incid 
%       Eather integer vector,  or structure with fields 'I' or 
%       'local' and 'imported'. 
%       Optionaly structure can have fields 'start_date' or 'dates'.
%       If 'dates' are given then start_date = dates(1).
%       'datesd' and 'start_date' must be  datetime class.
%
% Output:
%   res -- a structure with fields
%       res.local 
%       res.imported
%       res.I (= res.local  + res.imported)
%       res.start_date  (datetime object)
%       res.name
%       res.population

% 2020/07/28    Created

    if isstruct(incid)
        if isfield(incid,'local') && isfield(incid,'imported') ...
                && isfield(incid,'I')
            % probably incid created with this function
            res = incid;
            return
        elseif isfield(incid,'local') && isfield(incid,'imported') ...
                && ~isfield(incid,'I')
            validateattributes(incid.local, {'numeric'},...
                {'integer','vector','nonnegative'});
            validateattributes(incid.imported, {'numeric'},...
                {'integer','vector','nonnegative'});
            res.local = incid.local;
            res.imported = incid.imported;
        elseif isfield(incid,'I') && ...
                ~isfield(incid,'local') && ~isfield(incid,'imported')
            validateattributes(incid.I, {'numeric'},...
                {'integer','vector','nonnegative'});
            res.local = incid.I;
            res.local(1) = 0;
            res.imported = [incid.I(1),zeros(1,(length(incid.I)-1))];
        else
            error(['incid must be a vector or a structure with either',...
                'i)   a field called ''I'', or ',...
                'ii)  fields called ''local'' and ''imported'' or',...
                'iii) epiData object.']);
        end
        if isfield(incid,'start_date')
            res.start_date = datetime(incid.start_date);
            res.dates      = [];
        elseif isfield(incid,'dates')
            res.start_date = datetime(incid.dates(1));
        else
            res.start_date = [];
        end
        res.dates      = [];
    elseif isa(incid,'epiData')
        % epiData object, just copy data
        res.I = incid.I;
        res.local = incid.local;
        res.imported = incid.imported;
        res.start_date = incid.start_date;
        res.name  = incid.name;
        res.population = incid.population;
        return        
    else
        % vector
        validateattributes(incid, {'numeric'},...
            {'integer','vector','nonnegative'});
        res.local      = incid;
        res.imported   = zeros(1,length(incid));
        res.start_date = [];
        res.dates      = [];
    end

    if ~isrow(res.local)
        res.local = res.local';
    end
    if ~isrow(res.imported)
        res.imported = res.imported';
    end

    if res.local(1) > 0
        warning(['incid.local(1) is >0 but must be 0, as all cases on the first',...
            ' time step are assumed imported. This is corrected automatically',...
            ' by cases being transferred to incid.imported.'])
        res.imported(1) = res.imported(1) + res.local(1);
        res.local(1) = 0;
    end

    res.local(isnan(res.local)) = 0;
    res.imported(isnan(res.imported)) = 0;
    res.I = res.local + res.imported; % MB
    
    if isfield(incid,'name')
        res.name = incid.name;
    else
        res.name = [];
    end    
    if isfield(incid,'population')
        res.population = incid.population;
    else
        res.population = [];
    end

end

