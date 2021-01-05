%%% Estimation of the effective reproduction number over time (Rt) for 
function [Rt,Rt_adj,varargout]=wallinga_teunis(incid_data,serint_dist,varargin)

tmax=length(incid_data);

if length(serint_dist)>2
    gpar=gamfit(obsV); % If you have observed data on the serial interval distribution
else
gpar=serint_dist; % If these are the parameters of a gamma distribution
end
g=gampdf(1:tmax,gpar(1),gpar(2)); % Probability distribution for the serial interval

% Probability that a "true" case on day i was infected by another case on day j
p=zeros(tmax,tmax-1); % Initialize matrix to store values
for i=2:tmax % Loop through days 2 through tmax
    if incid_data(i)>0 % Only calculate p(i,j) if a case occurred on day i
        for j=1:i-1 % Loop through all days j BEFORE day i (i.e. from day 1 up to day i-1)
            if incid_data(j)>0 % Again, only calculate p(i,j) if a case occurred on day j
                p(i,j)=g(i-j)/(g((i-1):-1:1)*incid_data(1:(i-1))); % p(i,j) = probably the generation interval is (i-j) days over the sum of all possible generation intervals from cases occurring from day 1 of the epidemic up to day i-1
            end
        end
    end
end

% Calculate Rt for cases occurring during the first j = 1 to tmax-1 days of the epidemic 
% NOTE: Rt will be biased towards null for values approaching tmax  
Rt=zeros(tmax-1,1); % Initialize vector for R_t 
for j=1:tmax-1 % Loop through tmax-1 days
    Rt(j,1)=sum(p((j+1):end,j).*incid_data((j+1):end,1)); % R_t for day j is the sum of the probability that a case on day j infected a case on day j+1 UNTIL THE END OF THE EPIDEMIC times the number of cases ocurring on day j+1 until the end of the epidemic
end

% Adjustment for right-censoring of the data
Rt_adj=Rt./gamcdf((tmax-1:-1:1)',gpar(1),gpar(2));

% Calculate bootstrap CI
if ~isempty(varargin)
    nsamp=varargin{1};
    
    %Xt=zeros(tmax-1,nsamp);
    %for s=1:nsamp
    %    for j=1:tmax-1
    %        Xt(j,s)=sum(binornd(incid_data(j+1:end,1),p(j+1:end,j)));
    %    end
    %end
    
    Wt=gamcdf((tmax-1:-1:1)',gpar(1),gpar(2));
    %Xt_mean=mean(Xt,2)./Wt;
    %Xt_var=1./(Wt.^2).*var(Xt,0,2) + (1-Wt)./(Wt.^2).*mean(Xt,2);
    %Xt_lb=norminv(.025,Xt_mean,sqrt(Xt_var));
    %Xt_ub=norminv(.975,Xt_mean,sqrt(Xt_var));
    
    %Rt_adj=Xt_mean./incid_data(1:tmax-1,1);
    %Rt_ci=[Xt_lb./incid_data Xt_ub./incid_data(1:tmax-1,1)];
    
    
    %gmax=find(gamcdf(1:75,gpar(1),gpar(2))>.99,1); % Find the first day for which the CDF of the generation interval is >99%
    %incid_data=[incid_data; max(incid_data)*ones(gmax-1,1)]; % Pad the end of the incidence time series with values equal to the max number of cases previously observed on any given day. (This should be conservative if the epidemic is declining, but not necessarily if the epidemic is still increasing.)
    %g=gampdf(1:tmax+gmax,gpar(1),gpar(2)); % Probability distribution for the serial interval
    %p_any=zeros(tmax+gmax-1,tmax+gmax-2); % Initialize matrix to store values for the probability ANY case occurring on day j infected a case on day i
    p_any=zeros(tmax,tmax-1); % Initialize matrix to store values for the probability ANY case occurring on day j infected a case on day i
    for i=2:tmax % Loop through days 2 through tmax
        if incid_data(i,1)>0 % Only calculate p(i,j) if a case occurred on day i
            for j=1:i-1 % Loop through all days j BEFORE day i (i.e. from day 1 up to day i-1)
                if incid_data(j,1)>0 % Again, only calculate p(i,j) if a case occurred on day j
                    p_any(i,j)=(g(i-j)*incid_data(j,1))/(g((i-1):-1:1)*incid_data(1:(i-1),1)); % p(i,j) = probably the generation interval is (i-j) days over the sum of all possible generation intervals from ANY case occurring from day 1 of the epidemic up to day i-1
                end
            end
        end
    end

    %Y_samp=zeros(tmax+gmax-1,tmax+gmax-2);
    %Rt_samp=zeros(tmax+gmax-2,nsamp);
    Y_samp=zeros(tmax,tmax-1); 
    Yt_samp=zeros(nsamp,tmax-1);
    Rt_samp=zeros(tmax-1,nsamp);
    for s=1:nsamp
        for i=2:tmax
            if sum(p_any(i,:))>0
                Y_samp(i,:)=mnrnd(incid_data(i,1),p_any(i,:)); % Number of secondary cases on day i who were infected by any case occurring on days 1 to tmax
                %if i>tmax
                %    Y_samp(i,:)=mnrnd(round(rand*incid_data(i,1)),p_any(i,:));
                %end
            end
        end
        %Y_samp(tmax+1,:)=nbinrnd(sum(Y_samp),Wt);
        Yt_samp(s,:)=sum(Y_samp); %./Wt;
        %Rt_samp(:,s)=(sum(Y_samp)./Wt)'./incid_data(1:end-1,1);
    end
    %Rt_ci=prctile(Rt_samp(1:tmax-1,:)',[2.5 97.5])'; 
    Y_est=mean(Yt_samp)'./Wt;
    varYt=1./Wt.*(var(Yt_samp)')+(1-Wt)./(Wt.^2).*(mean(Yt_samp)');
    Y_est_ci=max(0,[norminv(.025,Y_est,sqrt(varYt)) norminv(.975,Y_est,sqrt(varYt))]);
    Rt_ci=Y_est_ci./(incid_data(1:end-1,1)*ones(1,2));
    varargout{1}=Rt_ci;
end

