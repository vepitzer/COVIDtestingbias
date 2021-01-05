function [T,T1] = exportRes(obj, varargin)
%EXPORT Export result to table and optionaly to .csv file
%
% Input Name-value pairs
%   exportRes(__,'file',name) 
%       save tables to csv file 'name'. Name should be without posfix .csv.
%       'name_R.csv' contain estimeted R, 'name_SI.csv' contains serial
%       interval.
%   exportRes(__,'delimiter',char)
%       set delimiter for csv file. Default is ';'.
%
% Output:
%   T --- table  
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

     
    if ~isempty(obj.IR)
        names = ["t_start", "t_end", "R", "IR","AR"];
        T = table(...
            obj.dates(obj.t_start)', ...
            obj.dates(obj.t_end)',...
            obj.mean_posterior', ...
            obj.IR',...
            obj.AR',....
            'VariableNames',names);
    else
        names = ["t_start", "t_end", "R"];
        T = table(...
            obj.dates(obj.t_start'), ...
            obj.dates(obj.t_end'),...
            obj.mean_posterior', ...
            'VariableNames',names);
    end
    
    if ~isempty(fnam)
        name = sprintf('%s_R.csv',fnam);
        writetable(T,name,'Delimiter',delimiter)
    end
    
end

