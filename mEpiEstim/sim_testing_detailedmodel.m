%%% PARAMETER VALUES %%%

tmax=110; %maximum number of time-steps for the simulation
dt=1; %time step = 1 day (for exponential distributed rates, approximate method)
dmax=round(tmax*dt); %maximum number of days for the simulation

beta1=[.85*ones(40,1); .3*ones(dmax-60,1); .5*ones(20,1)];  %transmission parameter
nu=1/2.9;  %rate of leaving latent period 
gamma=1/3.7;  %rate of recovery

%%% DEFINE TESTING SCENARIOS
nscen=8; % Number of scenarios

% Probability that a "true" symptomatic case is tested
Ptest_case_det=.05*ones(dmax,nscen); % Default value
Ptest_case_det(:,1)=[.05*ones(20,1); (.05:.0025:.15)'; .15*ones(dmax-61,1)]; % Small gradual increase in testing prob.
Ptest_case_det(:,2)=[.05*ones(30,1); .15*ones(dmax-30,1)]; % Small sudden increase in testing probability
Ptest_case_det(:,3)=[.05*ones(20,1); (.05:.0025:.25)'; .25*ones(dmax-101,1)]; % Large gradual increase in testing prob.
Ptest_case_det(:,4)=[.05*ones(30,1); .25*ones(dmax-30,1)]; % Large sudden increase in testing probability

% Probability that a severe case is tested
Ptest_severe=.5*ones(dmax,nscen); % Default value
Ptest_severe(:,2)=[.5*ones(30,1); (.5:.005:.75)'; .75*ones(dmax-81,1)]; % Also gradual increase in testing prob. for severe cases

% Ratio of the number tested to the number of "true" cases (up to day 50)
test_ratio_det=.2*ones(50,nscen); % Default value
test_ratio_det(:,[3 4])=.3*ones(50,2); % Higher testing capacity
test_ratio_det(:,5)=(.104:.004:.3)'; % Gradual increase in testing capacity
test_ratio_det(:,6)=(.3:-.004:.104)'; % Gradual decrease in testing capacity
test_ratio_det(:,7)=[.2*ones(30,1); .1*ones(20,1)]; % Abrupt decrease in testing capacity
test_ratio_det(:,8)=[.2*ones(30,1); .3*ones(20,1)]; % Abrupt increase in testing capacity

% Exponential growth rate of testing capacity (after day 50)
test_r=.016*ones(nscen,1); % Default value
test_r([4 7])=.02; % Faster growth in testing
test_r([3 8])=.013; % Slower growth in testing

%%% INITIALIZE OUTPUT VARIABLES %%%
NewCase_det=zeros(round(tmax*dt),nscen);
Symptomatic=zeros(round(tmax*dt)+25,nscen);
Severe=zeros(round(tmax*dt)+25,nscen);
Ntested_det=ones(round(tmax*dt)+25,nscen);
Nsymp=zeros(round(tmax*dt)+25,nscen);
Nsevere=zeros(round(tmax*dt)+25,nscen);
Npositive_det=zeros(round(tmax*dt)+25,nscen);
Ppositive_det=zeros(round(tmax*dt)+25,nscen);
Psevere=zeros(round(tmax*dt)+25,nscen);

R0_true_det=zeros(nscen,1);
R0_obs_det=zeros(nscen,1);
R0_tests_det=zeros(nscen,1);
R0ci_true_det=zeros(nscen,2);
R0ci_obs_det=zeros(nscen,2);
R0ci_tests_det=zeros(nscen,2);

Rt_true_det=zeros(round(tmax*dt)-1,nscen);
Rt_true_adj_det=zeros(round(tmax*dt)-1,nscen);
Rt_true_det_ci=zeros(round(tmax*dt)-1,2,nscen);
Rt_obs_det=zeros(round(tmax*dt)-1,nscen);
Rt_obs_adj_det=zeros(round(tmax*dt)-1,nscen);
Rt_obs_det_ci=zeros(round(tmax*dt)-1,2,nscen);
Rt_obs_full_det=zeros(round(tmax*dt)+24,nscen);
Rt_obs_full_adj_det=zeros(round(tmax*dt)+24,nscen);
Rt_tests_det=zeros(round(tmax*dt)-1,nscen);
Rt_tests_adj_det=zeros(round(tmax*dt)-1,nscen);
Rt_tests_det_ci=zeros(round(tmax*dt)-1,2,nscen);
Rt_tests_full_det=zeros(round(tmax*dt)+24,nscen);
Rt_tests_full_adj_det=zeros(round(tmax*dt)+24,nscen);

bias_upperci_det=zeros(dmax-1,nscen); 
bias_lowerci_det=zeros(dmax-1,nscen);

for scenario=1:nscen

%%% INITIALIZE STATE VARIABLES %%%
S=zeros(tmax,1);
E=zeros(tmax,1);
I=zeros(tmax,1);
R=zeros(tmax,1);
N=zeros(tmax,1);

Latent=zeros(tmax,1);
Infectious=zeros(tmax,1);
Recover=zeros(tmax,1);

N(1)=1000000;  %initial number of people 

E(1)=0;
I(1)=10;
R(1)=0;
S(1)=N(1)-E(1)-I(1)-R(1); 
Infectious(1)=I(1);

%%% START OF LOOP %%%
time=zeros(tmax,1);
for t=2:tmax
    
%%% RATES %%%
%%% Assuming exponentially distributed latent/infectious period
%%% (Approximate method)
Latent(t)=poissrnd(S(t-1)*beta1(ceil(t*dt))*I(t-1)/N(t-1)*dt); %Number of events is Poisson distributed if we assume they occur at a constant rate over short time step
Infectious(t)=poissrnd(nu*E(t-1)*dt);
Recover(t)=poissrnd(gamma*I(t-1)*dt);


%%% TRANSMISSION DYNAMICS %%%

%People
S(t)=max(S(t-1)-Latent(t),0);
E(t)=max(E(t-1)+Latent(t)-Infectious(t),0);
I(t)=max(I(t-1)+Infectious(t)-Recover(t),0);
R(t)=max(R(t-1)+Recover(t),0);
N(t)=S(t)+E(t)+I(t)+R(t);

%%% END OF LOOP
if I(t)==0, break, end

end
tf=t;

tstep=round(1/dt); %Number of time steps per day 
NewCase_det(1,scenario)=10; 
for i=2:round(tmax*dt)
    NewCase_det(i,scenario)=sum(Infectious(tstep*(i-1)+1:tstep*i,1)); %Number of new cases per day is the sum of individuals with symptom onset occurring during the 20 time steps that make up that day
end

%%% Incorporate testing practices and delays

symp_delay=gampdf(1:25,2.3,1); % from He et al
sev_delay=gampdf(1:25,5.1,.47); % Zhou et al

rep_dist=gampdf(1:25,1.85,3.57); % Default probability distribution for time from onset of infectiousness to testing/reporting 
rep_severe=gampdf(1:25,1,1); % reporting delay for severe cases

for i=1:round(tmax*dt)
    % Disease process
    for j=1:25
        Symptomatic(i+j,scenario)=Symptomatic(i+j,scenario)+binornd(binornd(NewCase_det(i,scenario),.6),symp_delay(j)); % Number of symptomatic cases is equal to the number of new infectious individuals on day i times the proportion symptomatic times the symptom delay distribution (pdf that symptoms appear j days after onset of infectiousness)
        Severe(i+j,scenario)=Severe(i+j,scenario)+binornd(binornd(Symptomatic(i,scenario),.2),sev_delay(j)); % Number of severe cases is equal to the number of symptomatic infections on day i times the proportion severe times the severe delay distribution (pdf of respiratory distress j days after onset of symptoms)
    end
    
    % Observation process
    Ntested_det(1,scenario)=1;
    if i<=50
    for j=1:25
        Ntested_det(i+j,scenario)=Ntested_det(i+j,scenario)+round(test_ratio_det(i,scenario)*NewCase_det(i,scenario)*rep_dist(j)); % Initial growth in number tested is equal to the number of new infections times the testing ratio times the reporting delay
    end 
    end
    Ntested_det(51:end,scenario)=round(Ntested_det(50,scenario)*exp(test_r(scenario)*(1:(dmax-25))')); % Subsequent growth in the number tested occurs exponentially with growth rate test_r
    for j=1:25
        Nsevere(i+j,scenario)=binornd(Severe(i,scenario),Ptest_severe(i,scenario)*rep_severe(j)); % Number of severe cases observed (test positive) with reporting delay
        Nsymp(i+j,scenario)=binornd(Symptomatic(i,scenario),Ptest_case_det(i,scenario)*rep_dist(j)); % Number of symptomatic cases observed (test positive) with reporting delay
        Npositive_det(i+j,scenario)=Npositive_det(i+j,scenario)+Nsevere(i+j,scenario); % Total number of positive tests is equal to the number positive plus the number of severe positive cases
        Npositive_det(i+j,scenario)=min(Ntested_det(i+j,scenario),Npositive_det(i+j,scenario)+Nsymp(i+j,scenario)); % Total number of positive tests is equal to the min of the number of tests or the number of severe + symptomatic positive cases
    end
end

Ppositive_det(:,scenario)=Npositive_det(:,scenario)./Ntested_det(:,scenario); % Proportion of tests positive
Psevere(:,scenario)=Nsevere(:,scenario)./(Nsymp(:,scenario)+Nsevere(:,scenario)); % Proportion of cases that are severe

% Estimate R0 from the growth rate of the epidemic (days 11-30)
[true_growthrate_det,~,true_r_stats_det]=glmfit(1:20,NewCase_det(21:40,scenario),'Poisson');
r_true_det=true_growthrate_det(2);
[obs_growthrate_det,~,obs_r_stats_det]=glmfit(1:20,Npositive_det(21:40,scenario),'Poisson');
r_obs_det=obs_growthrate_det(2);
[test_growthrate_det,~,test_r_stats_det]=glmfit(1:20,Ntested_det(21:40,scenario),'Poisson');
r_test_det=test_growthrate_det(2);

V=1/nu+1/gamma; % Duration of the serial interval (in days)
f=1/(V*nu); % Proportion of the serial interval in the latent period

R0_true_det(scenario,1)=r_true_det^2*f*(1-f)*V^2+r_true_det*V+1;
R0_obs_det(scenario,1)=r_obs_det^2*f*(1-f)*V^2+r_obs_det*V+1;
R0_tests_det(scenario,1)=r_test_det^2*f*(1-f)*V^2+r_test_det*V+1;

R0samp_true_det=zeros(1000,1);
R0samp_obs_det=zeros(1000,1);
R0samp_tests_det=zeros(1000,1);
for i=1:1000
    r=normrnd(r_true_det,true_r_stats_det.se(2));       
    R0samp_true_det(i,1)=r^2*f*(1-f)*V^2+r*V+1;
    r=normrnd(r_obs_det,obs_r_stats_det.se(2));       
    R0samp_obs_det(i,1)=r^2*f*(1-f)*V^2+r*V+1;
    r=normrnd(r_test_det,test_r_stats_det.se(2));       
    R0samp_tests_det(i,1)=r^2*f*(1-f)*V^2+r*V+1;
end
R0ci_true_det(scenario,:)=prctile(R0samp_true_det,[2.5 97.5]);
R0ci_obs_det(scenario,:)=prctile(R0samp_obs_det,[2.5 97.5]);
R0ci_tests_det(scenario,:)=prctile(R0samp_tests_det,[2.5 97.5]);

% Estimate Rt using the Wallinga-Teunis method
[Rt_true_det(:,scenario),Rt_true_adj_det(:,scenario),Rt_true_det_ci(:,:,scenario)]=wallinga_teunis(NewCase_det(1:round(tmax*dt),scenario),[4.79 1.357],100);
[Rt_obs_det(:,scenario),Rt_obs_adj_det(:,scenario),Rt_obs_det_ci(:,:,scenario)]=wallinga_teunis(Npositive_det(1:round(tmax*dt),scenario),[4.79 1.357],100);
[Rt_obs_full_det(:,scenario),Rt_obs_full_adj_det(:,scenario)]=wallinga_teunis(Npositive_det(:,scenario),[4.79 1.357]);
[Rt_tests_det(:,scenario),Rt_tests_adj_det(:,scenario),Rt_tests_det_ci(:,:,scenario)]=wallinga_teunis(Ntested_det(1:round(tmax*dt),scenario),[4.79 1.357],100);
[Rt_tests_full_det(:,scenario),Rt_tests_full_adj_det(:,scenario)]=wallinga_teunis(Ntested_det(:,scenario),[4.79 1.357]);

% Estimate the bias in estimates of Rt (as indicated by upper CI<0 or lower CI>0)
bias_upperci_det(:,scenario)=Rt_obs_det_ci(1:dmax-1,2,scenario)-Rt_true_adj_det(1:dmax-1,scenario);
bias_lowerci_det(:,scenario)=Rt_obs_det_ci(1:dmax-1,1,scenario)-Rt_true_adj_det(1:dmax-1,scenario);

% If bias is NaN, set equal to value from the previous day
for i=1:length(bias_upperci_det)
    if isnan(bias_upperci_det(i,scenario))
        if i==1
            bias_upperci_det(i,scenario)=0;
        else
            bias_upperci_det(i,scenario)=bias_upperci_det(i-1,scenario);
        end
    end
    if isnan(bias_lowerci_det(i,scenario))
        if i==1
            bias_lowerci_det(i,scenario)=0;
        else
            bias_lowerci_det(i,scenario)=bias_lowerci_det(i-1,scenario);
        end
    end
end    



% Estimate Rt using EpiEstim (method of Cori et al)
config = make_config('mean_si',6.5,'std_si',4); %'wind_len',3,

Rt_epiestim_true_det(scenario)=estimate_R(NewCase_det(:,scenario),'parametric_si',config);
Rt_epiestim_obs_det(scenario)=estimate_R(Npositive_det(:,scenario),'parametric_si',config);

Rt_true_det_epiestim(:,scenario)=Rt_epiestim_true_det(scenario).median_posterior';
Rt_true_det_epiestim_ci(:,:,scenario)=[Rt_epiestim_true_det(scenario).quantile_0025_posterior' Rt_epiestim_true_det(scenario).quantile_0975_posterior'];
Rt_obs_det_epiestim(:,scenario)=Rt_epiestim_obs_det(scenario).median_posterior';
Rt_obs_det_epiestim_ci(:,:,scenario)=[Rt_epiestim_obs_det(scenario).quantile_0025_posterior' Rt_epiestim_obs_det(scenario).quantile_0975_posterior'];

% Estimate the bias in estimates of Rt (as indicated by upper CI<0 or lower CI>0)
Rt_bias_det_epiestim(:,scenario)=max(0,Rt_obs_det_epiestim_ci(:,1,scenario)-Rt_true_det_epiestim(:,scenario));

end
%% PLOT RESULTS %%%

figure
for scenario=1:nscen
subplot(3,3,scenario)
hold on
bar(NewCase_det(:,scenario),'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
plot(Ntested_det(:,scenario),'r','LineWidth',2)
plot(Npositive_det(:,scenario),'b','LineWidth',2)
xlim([0 round(tmax*dt)+1])
ylim([0 15000])
xlabel('time (days)')
ylabel('# of cases')
title(join({'Scenario',num2str(scenario)}),'FontSize',11)
if scenario==nscen
    legend('"True" cases','Tests','Observed cases','Location','NW') 
end
end

%% Cases, percent positive, and Rt estimates

figure
for scenario=1:nscen
subplot(nscen,5,5*(scenario-1)+1)
%semilogy(cumsum(Npositive(:,scenario)),'b') 
%hold on
%semilogy(cumsum(Ntested(:,scenario)),'r')
%semilogy(cumsum(NewCase(:,scenario)),'--k')
hold on
[ax,y1,y2]=plotyy(1:dmax+25,Ntested_det(:,scenario)/100,1:dmax+25,Npositive_det(:,scenario)/100);
set(y1,'Color','r','LineWidth',1)
set(y2,'Color','b','LineWidth',1)
set(ax(1),'XLim',[20 dmax],'YLim',[0 100],'YTick',0:50:100,'YColor','k')
set(ax(2),'XLim',[20 dmax],'YLim',[0 15],'YTick',0:5:15,'YColor','k')
if scenario==1
legend('Test','Positive cases','Location','NW') 
end
ylabel('# tests/10K')

subplot(nscen,5,5*(scenario-1)+2)
hold on
plot(100*Ppositive_det(:,scenario),'Color',[.5 0 .5])
%plot(100*Psevere(:,scenario),'Color',[1 .5 0])
ylabel('% positive')
ylim([0 100])
xlim([20 dmax])
if scenario==1
    legend('Percent positive') %,'Percent severe')
end

subplot(nscen,5,5*(scenario-1)+3)
hold on
%plot(find(Rt_true_adj_det(:,scenario)),Rt_true_adj_det(Rt_true_adj_det(:,scenario)>0,scenario),'k')
%plot(find(Rt_obs_adj(:,scenario)),Rt_obs_adj(Rt_obs_adj(:,scenario)>0,scenario),'Color',[0 .7 0])
%fill([find(Rt_obs_adj(:,scenario)); flipud(find(Rt_obs_adj(:,scenario))); find(Rt_obs_adj(:,scenario),1)],[Rt_obs_det_ci(find(Rt_obs_adj(:,scenario)),1,scenario); Rt_obs_det_ci(flipud(find(Rt_obs_adj(:,scenario))),2,scenario); Rt_obs_det_ci(find(Rt_obs_adj(:,scenario),1),1,scenario)],[.9 1 .9],'EdgeColor',[.9 1 .9])
%fill([find(Rt_true_adj_det(:,scenario)); flipud(find(Rt_true_adj_det(:,scenario))); find(Rt_true_adj_det(:,scenario),1)],[Rt_true_det_ci(find(Rt_true_adj_det(:,scenario)),1,scenario); Rt_true_det_ci(flipud(find(Rt_true_adj_det(:,scenario))),2,scenario); Rt_true_det_ci(find(Rt_true_adj_det(:,scenario),1),1,scenario)],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
%plot(find(Rt_true_adj_det(:,scenario)),Rt_true_adj_det(Rt_true_adj_det(:,scenario)>0,scenario),'k')
%plot(find(Rt_obs_adj(:,scenario)),Rt_obs_adj(Rt_obs_adj(:,scenario)>0,scenario),'Color',[0 .7 0])
%plot(Rt_true_det(:,scenario),':k')
%plot(Rt_obs_det(:,scenario),':','Color',[0 .7 0])
plot(Rt_epiestim_true_det(1).t_end',Rt_true_det_epiestim(:,scenario),'k')
plot(Rt_epiestim_obs_det(1).t_end',Rt_obs_det_epiestim(:,scenario),'Color',[0 .7 0],'LineWidth',1)
fill([Rt_epiestim_obs_det(1).t_end'; flipud(Rt_epiestim_obs_det(1).t_end'); Rt_epiestim_obs_det(1).t_end(1)],[Rt_obs_det_epiestim_ci(:,1,scenario); flipud(Rt_obs_det_epiestim_ci(:,2,scenario)); Rt_obs_det_epiestim_ci(1,1,scenario)],[.9 1 .9],'EdgeColor',[.9 1 .9])
fill([Rt_epiestim_true_det(1).t_end'; flipud(Rt_epiestim_true_det(1).t_end'); Rt_epiestim_true_det(1).t_end(1)],[Rt_true_det_epiestim_ci(:,1,scenario); flipud(Rt_true_det_epiestim_ci(:,2,scenario)); Rt_true_det_epiestim_ci(1,1,scenario)],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
plot(Rt_epiestim_true_det(1).t_end',Rt_true_det_epiestim(:,scenario),'k')
plot(Rt_epiestim_obs_det(1).t_end',Rt_obs_det_epiestim(:,scenario),'Color',[0 .7 0],'LineWidth',1)
plot((0:dmax)',ones(dmax+1,1),'--k')
plot([20 20],[0 6],'k','LineWidth',.75)
ylabel('R_t')
if scenario==1
legend('True','Observed'); 
end
xlim([20 70])
ylim([0 8])
if scenario==8
    xlabel('Days since start of epidemic')
end

subplot(nscen,5,5*(scenario-1)+4)
hold on
%plot(find(Rt_true_adj_det(:,scenario)),Rt_true_adj_det(Rt_true_adj_det(:,scenario)>0,scenario),'k')
%plot(find(Rt_obs_adj(:,scenario)),Rt_obs_adj(Rt_obs_adj(:,scenario)>0,scenario),'Color',[0 .7 0])
%fill([find(Rt_obs_adj(:,scenario)); flipud(find(Rt_obs_adj(:,scenario))); find(Rt_obs_adj(:,scenario),1)],[Rt_obs_det_ci(find(Rt_obs_adj(:,scenario)),1,scenario); Rt_obs_det_ci(flipud(find(Rt_obs_adj(:,scenario))),2,scenario); Rt_obs_det_ci(find(Rt_obs_adj(:,scenario),1),1,scenario)],[.9 1 .9],'EdgeColor',[.9 1 .9])
%fill([find(Rt_true_adj_det(:,scenario)); flipud(find(Rt_true_adj_det(:,scenario))); find(Rt_true_adj_det(:,scenario),1)],[Rt_true_det_ci(find(Rt_true_adj_det(:,scenario)),1,scenario); Rt_true_det_ci(flipud(find(Rt_true_adj_det(:,scenario))),2,scenario); Rt_true_det_ci(find(Rt_true_adj_det(:,scenario),1),1,scenario)],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
%plot(find(Rt_true_adj_det(:,scenario)),Rt_true_adj_det(Rt_true_adj_det(:,scenario)>0,scenario),'k')
%plot(find(Rt_obs_adj(:,scenario)),Rt_obs_adj(Rt_obs_adj(:,scenario)>0,scenario),'Color',[0 .7 0])
%plot(Rt_true_det(:,scenario),':k')
%plot(Rt_obs_det(:,scenario),':','Color',[0 .7 0])
plot(Rt_epiestim_true_det(1).t_end',Rt_true_det_epiestim(:,scenario),'k')
plot(Rt_epiestim_obs_det(1).t_end',Rt_obs_det_epiestim(:,scenario),'Color',[0 .7 0],'LineWidth',1)
fill([Rt_epiestim_obs_det(1).t_end'; flipud(Rt_epiestim_obs_det(1).t_end'); Rt_epiestim_obs_det(1).t_end(1)],[Rt_obs_det_epiestim_ci(:,1,scenario); flipud(Rt_obs_det_epiestim_ci(:,2,scenario)); Rt_obs_det_epiestim_ci(1,1,scenario)],[.9 1 .9],'EdgeColor',[.9 1 .9])
fill([Rt_epiestim_true_det(1).t_end'; flipud(Rt_epiestim_true_det(1).t_end'); Rt_epiestim_true_det(1).t_end(1)],[Rt_true_det_epiestim_ci(:,1,scenario); flipud(Rt_true_det_epiestim_ci(:,2,scenario)); Rt_true_det_epiestim_ci(1,1,scenario)],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
plot(Rt_epiestim_true_det(1).t_end',Rt_true_det_epiestim(:,scenario),'k')
plot(Rt_epiestim_obs_det(1).t_end',Rt_obs_det_epiestim(:,scenario),'Color',[0 .7 0],'LineWidth',1)
plot((0:dmax)',ones(dmax+1,1),'--k')
plot([20 20],[0 6],'k','LineWidth',.75)
%ylabel('R_t')
if scenario==1
legend('True','Observed') 
end
xlim([70 dmax-4])
ylim([0 2])

subplot(nscen,5,5*scenario)
hold on
bar(find(bias_lowerci_det(:,scenario)>0),bias_lowerci_det(bias_lowerci_det(:,scenario)>0,scenario),'k')
bar(find(bias_upperci_det(:,scenario)<0),bias_upperci_det(bias_upperci_det(:,scenario)<0,scenario),'k')
if scenario==1
legend('Bias')
end
xlim([20 dmax-1])
ylim([-2 2])

end

%%
gtext('A)','FontWeight','bold')
gtext('B)','FontWeight','bold')
gtext('C)','FontWeight','bold')
gtext('D)','FontWeight','bold')
gtext('E)','FontWeight','bold')
gtext('F)','FontWeight','bold')
gtext('G)','FontWeight','bold')
gtext('H)','FontWeight','bold')

%% Estimated R0 values
figure
hold on
plot(.9:1:(nscen-.1),R0_true_det,'ok')
plot(1:nscen,R0_obs_det,'ob')
plot(1.1:(nscen+.1),R0_tests_det,'or')
plot([.9:1:(nscen-.1); .9:1:(nscen-.1)],R0ci_true_det','k')
plot([1:nscen; 1:nscen],R0ci_obs_det','b')
plot([1.1:1:(nscen+.1); 1.1:1:(nscen+.1)],R0ci_tests_det','r')
plot([0 nscen+1],[mean(R0_true_det) mean(R0_true_det)],'Color',[.7 .7 .7])
legend('True value','Confirmed cases','Number of tests','95% CI')
ylabel('R_0')
xlabel('Scenario')
set(gca,'XTick',1:nscen)
