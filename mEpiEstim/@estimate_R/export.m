function [T,T1] = export(obj, varargin)
%EXPORT Export result to table and optionaly to .csv file
%
% Input Name-value pairs
%   export__,'file',name) 
%       save tables to csv file 'name'. Name should be without posfix .csv.
%       'name_R.csv' contain estimeted R, 'name_SI.csv' contains serial
%       interval.
%   export(__,'delimiter',char)
%       set delimiter for csv file. Default is ';'.
%
% Output:
%   T --- estimated R
%   T1 -- SI
%

    % Default values
    fnam = [];
    delimiter = ';';
    
    % Scan input
    if ~isempty(varargin)
        for n = 1:2:length(varargin)
            if isempty(varargin{n+1})
                % skip blanks
                continue
            end
            switch lower(varargin{n})
                case {'delimiter'}
                    delimiter = varargin{n + 1};                
                case {'file','save'}
                    fnam = varargin{n + 1};
                otherwise
                    error('Unknown property %s.',varargin{n})
            end
        end
    end

    % Estimated R
    names = ["Var1",...
        "t_start", "t_end", "Mean(R)", "Std(R)",...
        "Quantile.0.025(R)", "Quantile.0.05(R)", "Quantile.0.25(R)",...
        "Median(R)", "Quantile.0.75(R)", "Quantile.0.95(R)",...
        "Quantile.0.975(R)"];
     
    T = table(...
        num2str((1:length(obj.config.t_start))'),...
        obj.config.t_start', ...
        obj.config.t_end',...
        obj.mean_posterior', ...
        obj.std_posterior',...
        obj.quantile_0025_posterior',...
        obj.quantile_005_posterior',...
        obj.quantile_025_posterior',...
        obj.median_posterior',...
        obj.quantile_075_posterior',...
        obj.quantile_095_posterior',...
        obj.quantile_0975_posterior',...
        'VariableNames',names);
    
    if ~isempty(fnam)
        name = sprintf('%s_R.csv',fnam);
        writetable(T,name,'Delimiter',delimiter)
    end
    
    % Serial interval    
    if isvector(obj.config.si_distr)
        names1 = [" ","x"];
        % Default values
        
        T1 = table(...
            num2str((1:length(obj.config.si_distr))'),...
            obj.config.si_distr',...
            'VariableNames',names1);
        
        if ~isempty(fnam)
            name = sprintf('%s_SI.csv',fnam);
            writetable(T1,name,'Delimiter',delimiter)
        end
    else
        %TODO
        T1 = [];
    end
end

