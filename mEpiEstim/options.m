classdef options < handle
    %OPTIONS Create empty object for estim_R plot options

    properties
        col                  % color
        transp               % transparentcity
        xlim                 % x-axis limits
        ylim                 % y-axis limits
        xlab                 % x-axis label
        ylab                 % y-axis label
        title                % graph title
        extend_title  = true % include additional data in title ?
        plot_type            % set plot type: line, bar, stairs, scatter,...
        scale                % linear, semilogy
    end
    
    % Incidence properties
    properties
        interval = 1
    end
    
    % SI
    properties
        prob_min = 0.001
    end
    
    methods
        function obj = options(varargin)
            %OPTIONS Construct an instance of this class
            
            if isempty(varargin)
                return
            end
            
            for n = 1:2:length(varargin)
                if isempty(varargin{n+1})
                    continue
                end
                switch lower(varargin{n})
                    case 'col'
                        obj.col = varargin{n+1};                    
                    case 'extend_title'
                        obj.extend_title = chkOnOff(varargin{n+1});
                    case 'plot_type'
                        obj.plot_type = lower(varargin{n+1});
                    case 'scale'
                        obj.scale = lower(varargin{n+1});                        
                    case 'title'
                        obj.title = varargin{n+1};                        
                    case 'xlab'
                        obj.xlab = varargin{n+1};
                    case 'xlim'
                        obj.xlim = varargin{n+1};    
                    case 'ylab'
                        obj.ylab = varargin{n+1};
                    case 'ylim'
                        obj.ylim = varargin{n+1};                         
                end
            end
        end
        function set.plot_type(obj,plot_type)
            switch plot_type
                case {'bar','line','scatter','stairs','area'}
                    obj.plot_type = plot_type;
                otherwise
                    error('Invalid plot type %s',plot_type)
            end
            
        end
        function set.prob_min(obj,prob_min)
            validateattributes(prob_min, {'numeric'},...
                {'real','scalar','positive','<',1});
            obj.prob_min = prob_min;
        end        
        function set.title(obj,title)
            validateattributes(title, {'char'},{'vector'});
            obj.title = title;
        end         
        function set.transp(obj,transp)
            validateattributes(transp, {'numeric'},...
                {'real','scalar','>=',0,'<=',1});
            obj.transp = transp;
        end         
        function set.xlim(obj,xlim)
            validateattributes(xlim, {'numeric'},...
                {'real','vector','increase','numel',2});
            obj.xlim = xlim;
        end
        function set.ylim(obj,ylim)
            validateattributes(ylim, {'numeric'},...
                {'real','vector','increasing','numel',2});
            obj.ylim = ylim;
        end        
    end
end

