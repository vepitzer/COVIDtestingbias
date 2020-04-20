%%% Estimation of the effective reproduction number over time (Rt) for 
function [Rt,Rt_adj]=wallinga_teunis(incid_data,serint_dist)

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
Rt=zeros(tmax,1); % Initialize vector for R_t 
for j=1:tmax-1 % Loop through tmax-1 days
    Rt(j,1)=sum(p((j+1):end,j).*incid_data((j+1):end,1)); % R_t for day j is the sum of the probability that a case on day j infected a case on day j+1 UNTIL THE END OF THE EPIDEMIC times the number of cases ocurring on day j+1 until the end of the epidemic
end

% Adjustment for right-censoring of the data
Rt_adj=Rt./gamcdf((tmax:-1:1)',gpar(1),gpar(2));


