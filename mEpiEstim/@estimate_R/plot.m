function plot(x,varargin)
%PLOT_ESTIMATE_R Plot outputs of estimate_r

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on R function plot from EpiEstim pacgage
% https://cran.r-project.org/web/packages/EpiEstim/index.html
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2020/07/29    Created

% Input:
%   x  
%       The output of function estimate_R

% Optional input:
%   what
%       plot(__,'what',{'all','incid','R','SI','IR'})
%       plot(__,'what',['incid','R',...]) -- plots incidence and R  etc
%   add_imported_cases
%       plot(__,'add_imported_cases',true) -- add the incidence of imported cases.
%
%   options_I
%
%   options_R 
%
%   options_SI
%
%   options_AR
%
%   legend
%       plot(__,'legend',true) -- add legends on the plots (true by default)
%   grid
%       plot(__,'grid',true) -- add grid on the plots (true by default)
%   title
%       plot(__,'title',str) -- add graph title
% 
%

    lang = '';
    
    % Default values
    what = [];
    add_imported_cases = false;
    legendf = true;
    gridf   = true;
    titlStr = [];
    
    options_I = [];
    options_R = [];
    options_SI = [];
    options_IR = [];
    
    % scan options
    if ~isempty(varargin)
        for n = 1:2:length(varargin)
            if isempty(varargin{n+1})
                % skip blanks
                continue
            end
            switch lower(varargin{n})
                case 'legend'
                    legendf = chkOnOff(varargin{n+1});
                case 'grid'
                    gridf = chkOnOff(varargin{n+1});    
                case 'title'
                    titlStr = varargin{n+1};
                case 'what'
                    % what to plot
                    tmp = lower(varargin{n+1});
                    if strcmpi(tmp,'all')
                        what = [1,1,1,0];
                    else
                         what = zeros(1,4);
                        if contains(tmp,'r')
                            what(2) = 1;
                        end
                        if contains(tmp,'incid')
                            what(1) = 1;
                        end
                        if contains(tmp,'si')
                            what(3) = 1;
                        end
                        if contains(tmp,'ir')
                            what(4) = 1;
                        end                        
                    end
                case 'add_imported_cases'
                    add_imported_cases = varargin{n+1};
                case {'opt_r','options_r'}
                    options_R = varargin{n+1};
                    if ~isa(options_R,'options')
                        error('Exepted input must be ''options'' object.')
                    end                             
                case {'opt_i','options_i'}
                    options_I = varargin{n+1};
                    if ~isa(options_I,'options')
                        error('Exepted input must be ''options'' object.')
                    end                      
                case {'opt_si','options_si'}
                    options_SI = varargin{n+1};
                    if ~isa(options_SI,'options')
                        error('Exepted input must be ''options'' object.')
                    end        
                case {'opt_ir','options_ir'}
                    options_IR = varargin{n+1};
                    if ~isa(options_IR,'options')
                        error('Exepted input must be ''options'' object.')
                    end                      
            end
        end
    end

    % check
    if isempty(what)
        what = [1,1,1,0];
    end
    
    if isempty(options_I)
        options_I = options();
    end
    if isempty(options_I.transp)
        options_I.transp = 0.7;
    end
    if isempty(options_I.xlab)
        options_I.xlab = 'Time';
    end
    if isempty(options_I.ylab)
        options_I.ylab = 'Incidence';
    end
    if isempty(options_I.title)
        options_I.title = 'Epidemic curve';
    end
    if isempty(options_I.plot_type)
        options_I.plot_type   = 'bar';
    end
    
    if isempty(options_R)
        options_R = options();
    end
    if isempty(options_R.col)
        options_R.col = 0.8*[1,1,1];
    end
    if isempty(options_R.transp)
        options_R.transp = 0.2;
    end
    if isempty(options_R.xlab)
        options_R.xlab = 'Time';
    end
    if isempty(options_R.ylab)
        options_R.ylab = 'R';
    end
    if isempty(options_R.title)
        options_R.title = 'Reproduction Number';
    end
    
    if isempty(options_SI)
        options_SI = options();
    end
    if isempty(options_SI.col)
        options_SI.col = 'black';
    end
    if isempty(options_SI.transp)
        options_SI.transp = 0.25;
    end
    if isempty(options_SI.xlab)
        options_SI.xlab = 'Time';
    end
    if isempty(options_SI.ylab)
        options_SI.ylab = 'Frequency';
    end
    if isempty(options_SI.title)
        options_SI.title = 'Explored SI distribution';
    end
    if isempty(options_SI.plot_type)
        options_SI.plot_type   = 'line';
    end
    
    if isempty(options_IR)
        options_IR = options();
    end
    if isempty(options_IR.col)
        options_IR.col = 'black';
    end
    if isempty(options_IR.transp)
        options_IR.transp = 0.25;
    end
    if isempty(options_IR.xlab)
        options_IR.xlab = 'Time';
    end
    if isempty(options_IR.ylab)
        options_IR.ylab = 'IR';
    end
    if isempty(options_IR.title)
        options_IR.title = sprintf('Incidence Rate per %d persons',x.N);
    end
    if isempty(options_IR.plot_type)
        options_IR.plot_type   = 'line';
    end
    
    % Extract data
    t_start = x.t_start;
    t_end   = x.t_end;
    mean_posterior          = x.mean_posterior;
    quantile_0025_posterior = x.quantile_0025_posterior;
    quantile_0975_posterior = x.quantile_0975_posterior;
   % method   = x.method;
    si_distr = x.si_distr;
    incid.I        = x.I;
    incid.local    = x.I_local;
    incid.imported = x.I_imported;
    
    % Limits
    T = length(incid.local);
    if ~isempty(x.dates)
        dates = x.dates;
        if isdatetime(dates)
            dmin  = datetime(datestr(min(datenum(x.dates) - 1)));
            dmax  = datetime(datestr(max(datenum(x.dates) + 1)));
        else
            dmin = 0;
            dmax = T + 1;
        end
    else
        dates = 1:T;
        dmin  = 0;
        dmax  = T + 1;
    end
  
    % start plot

    figure
    nplt = sum(what);
    switch nplt
        case 1
        case 2
            set(gcf,'Position',[50 50 832 624])
        case 3
            set(gcf,'Position',[50 50 832 832])
        case 4
            set(gcf,'Position',[50 50 832 832])            
    end
    np = 0;
    if isempty(titlStr)
        % no extra title
        tflag = true;
    else
        tflag = false;
    end
    
    %=====================================================================%
    % shows the epidemic curve
    %=====================================================================%
    if what(1) == 1
        np = np + 1;
        subplot(nplt,1,np)
        hold on
        if add_imported_cases
            switch options_I.plot_type
                case 'bar'
                    bar(dates,[incid.local;incid.imported],'stacked')
                case 'scatter'
                    scatter(dates,incid.local,[],'b')
                    scatter(dates,incid.imported,[],'r')                    
                case 'stairs'
                    stairs(dates,incid.local,'b','LineWidth',2)
                    stairs(dates,incid.imported,'r','LineWidth',2)
                otherwise
                    plot(dates,incid.local,'b','LineWidth',2)
                    plot(dates,incid.imported,'r','LineWidth',2)
            end
            legend('local','imported','Location','best')
        else
            switch options_I.plot_type
                case 'bar'
                    bar(dates,incid.I)
                case 'scatter'
                    scatter(dates,incid.I,[],'k')                
                case 'stairs'
                    stairs(dates,incid.I,'k','LineWidth',2)  
                otherwise
                    plot(dates,incid.I,'k','LineWidth',2)                    
            end
        end
        xlabel(options_I.xlab)
        ylabel(options_I.ylab)
        if ~isempty(options_I.xlim)
            xlim(options_I.xlim)
        else
            xlim([dmin,dmax])
        end
        if ~isempty(options_I.ylim)
            ylim(options_I.ylim)
        end
        if ~tflag
            txt = sprintf('%s\n%s',titlStr,options_I.title);
            tflag = true;
        else
            txt = options_I.title;
        end
        title(txt)
        if gridf
            grid on
        end
        hold off
    end
    
    %=====================================================================%
    % The first one shows the epidemic curve. The second one
    %shows the posterior mean and 95\% credible interval of the reproduction
    %number. The estimate for a time window is plotted at the end of the time
    %window.
    %=====================================================================%
    if what(2) == 1
        np = np + 1;
        subplot(nplt,1,np)
        hold on
        
        time_points = [];
        n = 0;
        while true
            n = n + 1;
            time_points(n) = length(t_start(n):t_end(n) - 1);
            if n == length(t_start)
                break
            end
        end
        
        if length(time_points) == length(unique(time_points))
            sliding_windows = false;
        else
            sliding_windows = true;
        end
        
        if sliding_windows
            t  = dates(t_end);
            t2 = [t,fliplr(t)];
            fa = [quantile_0025_posterior,fliplr(quantile_0975_posterior)];
            fill(t2,fa,options_R.col,'EdgeColor','none')
            plot(t,mean_posterior,'k','LineWidth',2)
        else
            if isdatetime(dates)
                t1 = datetime(datestr(datenum(dates(t_start) - 0.5)));
                t2 = datetime(datestr(datenum(dates(t_end) + 0.5)));
            else
                t1 = t_start - 0.5;
                t2 = t_end + 0.5;
            end
            
           % cola = colormap(jet(length(t_start)));
            for n = 1:length(t_start)
                t = [t1(n) t2(n) t2(n) t1(n)];
                fa = [quantile_0025_posterior(n) quantile_0025_posterior(n) ...
                    quantile_0975_posterior(n) quantile_0975_posterior(n)];
                %fill(t,fa,cola(n,:),'EdgeColor','none')
                 fill(t,fa,options_R.col,'EdgeColor','none')
                plot([t1(n) t2(n)],[mean_posterior(n) mean_posterior(n)],...
                    'k','LineWidth',2)
            end
        end
        if legendf
            legend('95%CrI','Mean','Location','best')
        end
        h = plot([dmin,dmax],[1,1],'r');
        h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        xlabel(options_R.xlab)
        ylabel(options_R.ylab)
        if ~isempty(options_R.xlim)
            xlim(options_R.xlim)
        else
            xlim([dmin,dmax])
        end
        if ~isempty(options_R.ylim)
            ylim(options_R.ylim)
        end
        if ~options_R.extend_title
            txt = options_R.title;
        else
            if isdatetime(x.dates)
                txt = sprintf('%s (%s R = %g)', options_R.title,...
                    x.dates(end),round(mean_posterior(end),2));
            else
                txt = sprintf('%s (R = %g)', options_R.title,...
                    round(mean_posterior(end),2));
            end
        end
        if ~tflag
            txt = sprintf('%s\n%s',titleStr,txt);
            tflag = true;
        end
        title(txt)
        ylabel(options_R.ylab)
        if gridf
            grid on
        end
        hold off
    end
    
    %=====================================================================%
    % shows the discrete distribution(s) of the serial interval.
    %=====================================================================%
    if what(3) == 1
        np = np + 1;
        subplot(nplt,1,np)        
        hold on
        if isvector(si_distr)
            tmp = cumsum(si_distr(si_distr >= options_SI.prob_min));
            si_distr_for_plot = si_distr(1:length(tmp));
            yymax = max(si_distr_for_plot);
            tsi = 0:length(si_distr_for_plot) - 1;
            switch options_SI.plot_type
                case 'bar'
                    bar( tsi, si_distr_for_plot)
                case 'scatter'
                    scatter(tsi,si_distr_for_plot,[],'k')
                case 'stairs'
                    stairs(tsi,si_distr_for_plot,'k','LineWidth',2)
                otherwise
                    plot( tsi, si_distr_for_plot,...
                        'k','LineWidth',2)
            end
        else
            %TODO
            switch options_SI.plot_type
                case 'area'
                    mean_R = mean(si_distr,1);
                    q005_R = quantile(si_distr,0.05,1);
                    q095_R = quantile(si_distr,0.95,1);
                    tmp = cumsum(mean_R(mean_R >= options_SI.prob_min));
                    mean_R = mean_R(1:length(tmp));
                    q005_R = q005_R(1:length(tmp));
                    q095_R = q095_R(1:length(tmp));
                    xxlim = length(tmp);
                    yymax = max(q095_R);
                    t = 0:length(mean_R)-1;
                    t2 = [t,fliplr(t)];
                    fa = [q005_R,fliplr(q095_R)];
                    fill(t2,fa,options_R.col,'EdgeColor','none')
                    plot(t,mean_R,'k','LineWidth',2)
                    if legendf
                        legend('95%CrI','Mean','Location','best')
                    end
                otherwise
                    xxlim = 0;
                    yymax = 0;
                    for k = 1:size(si_distr,1)
                        tmp = cumsum(si_distr(k,si_distr(k,:) >= options_SI.prob_min));
                        si_distr_for_plot = si_distr(k,1:length(tmp));
                        xxlim = max([xxlim,length(si_distr_for_plot)]);
                        yymax = max([yymax,max(si_distr_for_plot)]);
                        tsi = 0:length(si_distr_for_plot) - 1;
                        plot( tsi, si_distr_for_plot,...
                            'g','LineWidth',2)
                    end
%                     mean_R = mean(si_distr,1);
%                     q005_R = quantile(si_distr,0.05,1);
%                     q095_R = quantile(si_distr,0.95,1);
%                     tsi = 0:length(mean_R)-1;
%                     plot(tsi,mean_R,'k','LineWidth',2)
%                     plot(tsi,q005_R,'r','LineWidth',2)
%                     plot(tsi,q095_R,'r','LineWidth',2)
            end
        end
        xlabel(options_SI.xlab)
        ylabel(options_SI.ylab)
        if ~isempty(options_SI.xlim)
            xlim(options_SI.xlim)
        else
            if isvector(si_distr)
                xlim([0,length(si_distr_for_plot)-1]);
            else
                xlim([0,xxlim]);
            end
        end
        if ~isempty(options_SI.ylim)
            ylim(options_SI.ylim)
        end
        % plot SD
        if ~isempty(x.config.mean_si)
            xx = x.config.mean_si;
            h = plot( xx*[1,1],[0,yymax],'r');
            h.Annotation.LegendInformation.IconDisplayStyle = 'off';
        end
        if ~options_SI.extend_title
            txt = options_SI.title;
        else
            if ~isempty(x.config.mean_si)
                txt = sprintf('%s (mean = %g SD = %g)', options_SI.title, ...
                    round(x.config.mean_si,2),round(x.config.std_si,2));
            else
                txt = options_SI.title;
            end
        end
        if ~tflag
            txt = sprintf('%s\n%s',titleStr,txt);
            tflag = true;
        end
        title(txt)      
        if gridf
            grid on
        end
        hold off
    end
    
    %=====================================================================%
    % shows the incidence rate
    %=====================================================================%
    if what(4) == 1 && ~isempty(x.IR)
        np = np + 1;
        subplot(nplt,1,np)  
        hold on
%         if unique(t_end - t_start + 1) == 14
%             wnd = 'two-weeks';
%         elseif unique(t_end - t_start + 1) == 7
%             wnd = 'week';
%         else
            nwnd = unique(t_end - t_start + 1); 
            if isscalar(nwnd)
                wnd = sprintf('%d days',nwnd);
            else
                wnd = [];
            end
%         end
        IR = x.IR; 
        switch options_IR.plot_type
            case 'bar'
                bar(dates(t_end),IR)
            case 'scatter'
                scatter(dates(t_end),IR,[],'k')
            case 'stairs'
                stairs(dates(t_end),IR,'k','LineWidth',2)
            otherwise
                plot(dates(t_end),IR,'k','LineWidth',2)
        end
        if nwnd == 14
            if max(IR) > 10
                plot([dates(1),dates(end)],10*[1,1],'g','LineWidth',2)
            end
            if max(IR) > 40
                plot([dates(1),dates(end)],40*[1,1],'r','LineWidth',2)
            end
        end
        xlabel(options_IR.xlab)
        ylabel(options_IR.ylab)
        if ~isempty(wnd)
            txt = sprintf('%s per %s',options_IR.title,wnd);
        else
            txt = options_IR.title;
        end
        if options_IR.extend_title
            if isdatetime(x.dates)
                txt = sprintf('%s (%s IR = %g)', txt,...
                    x.dates(end),round(IR(end),0));
            else
                txt = sprintf('%s (IR = %g)', txt,...
                    round(IR(end),0));
            end
        end
        if ~tflag
            txt = sprintf('%s\n%s',titleStr,txt);
            %tflag = true;
        end
        title(txt)
        if ~isempty(options_IR.xlim)
            xlim(options_IR.xlim)
        else
            xlim([dmin,dmax])
        end
        if ~isempty(options_IR.ylim)
            ylim(options_IR.ylim)
        end
        if gridf
            grid on
        end
    end
    
end

