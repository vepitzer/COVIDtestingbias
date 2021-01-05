%%% PARAMETER VALUES %%%

tmax=1400; %maximum number of time-steps for the simulation
dt=0.05; %time step = 1/20 day (for exponential distributed rates, approximate method)

beta1=[.85*ones(round(tmax*dt)-20,1); .5*ones(20,1)];  %transmission parameter
nu=1/2.9;  %rate of leaving latent period 
gamma=1/3.7;  %rate of recovery

%%% DEFINE TESTING SCENARIOS

rep_dist_pars=gampdf(1:25,1.85,3.57); % Default probability distribution for time from onset of infectiousness to testing/reporting 

Ptest_case=zeros(round(tmax*dt),11);
Ptest_case(:,1)=.1*ones(round(tmax*dt),1); % Probability that a "true" case is tested
Ptest_case(:,2)=[.05*ones(round(tmax*dt)-51,1); (.05:.005:.25)'; .25*ones(10,1)];
Ptest_case(:,3)=[.25*ones(round(tmax*dt)-51,1); (.25:-.005:.05)'; .05*ones(10,1)];
Ptest_case(:,4:5)=.1*ones(round(tmax*dt),2); % Probability that a "true" case is tested
Ptest_case(:,6)=[.05*ones(round(tmax*dt)-30,1); .25*ones(30,1)];
Ptest_case(:,7)=[.25*ones(round(tmax*dt)-30,1); .05*ones(30,1)];
Ptest_case(:,8:11)=.1*ones(round(tmax*dt),4); % Probability that a "true" case is tested

test_ratio=zeros(round(tmax*dt),11);
test_ratio(:,1:3)=.5*ones(round(tmax*dt),3); % Ratio of the number tested to the number of "true" cases
test_ratio(:,4)=[.2*ones(round(tmax*dt)-51,1); (.2:.015:.8)'; .8*ones(10,1)]; % Ratio of the number tested to the number of "true" cases
test_ratio(:,5)=[.8*ones(round(tmax*dt)-51,1); (.8:-.015:.2)'; .2*ones(10,1)]; % Ratio of the number tested to the number of "true" cases
test_ratio(:,6:7)=.5*ones(round(tmax*dt),2); % Ratio of the number tested to the number of "true" cases
test_ratio(:,8)=[.2*ones(round(tmax*dt)-30,1); .8*ones(30,1)]; % Ratio of the number tested to the number of "true" cases
test_ratio(:,9)=[.8*ones(round(tmax*dt)-30,1); .2*ones(30,1)]; % Ratio of the number tested to the number of "true" cases
test_ratio(:,10:11)=.5*ones(round(tmax*dt),2); % Ratio of the number tested to the number of "true" cases

%%% INITIALIZE OUTPUT VARIABLES %%%
NewCase=zeros(round(tmax*dt),11);
Ntested=ones(round(tmax*dt)+25,11);
Npositive=zeros(round(tmax*dt)+25,11);
Ppositive=zeros(round(tmax*dt)+25,11);

R0_true=zeros(11,1);
R0_obs=zeros(11,1);
R0_tests=zeros(11,1);
R0ci_true=zeros(11,2);
R0ci_obs=zeros(11,2);
R0ci_tests=zeros(11,2);
Rt_true=zeros(round(tmax*dt)-1,11);
Rt_true_adj=zeros(round(tmax*dt)-1,11);
Rt_true_ci=zeros(round(tmax*dt)-1,2,11);
Rt_obs=zeros(round(tmax*dt)-1,11);
Rt_obs_adj=zeros(round(tmax*dt)-1,11);
Rt_obs_ci=zeros(round(tmax*dt)-1,2,11);
Rt_obs_full=zeros(round(tmax*dt)+24,11);
Rt_obs_full_adj=zeros(round(tmax*dt)+24,11);
Rt_tests=zeros(round(tmax*dt)-1,11);
Rt_tests_adj=zeros(round(tmax*dt)-1,11);
Rt_tests_ci=zeros(round(tmax*dt)-1,2,11);
Rt_tests_full=zeros(round(tmax*dt)+24,11);
Rt_tests_full_adj=zeros(round(tmax*dt)+24,11);

bias_upperci=zeros(round(tmax*dt)-1,11);
bias_lowerci=zeros(round(tmax*dt)-1,11);

for scenario=1:11

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

tstep=round(1/dt); %Number of time steps per day (= 1/.05 = 20)
NewCase(1,scenario)=10; 
for i=2:round(tmax*dt)
    NewCase(i,scenario)=sum(Infectious(tstep*(i-1)+1:tstep*i,1)); %Number of new cases per day is the sum of individuals with symptom onset occurring during the 20 time steps that make up that day
end


%%% Incorporate testing practices and delays

for i=1:round(tmax*dt)
    if scenario==10
    if i<30
        rep_dist=gampdf(1:25,1.85,3.57); % Probability distribution for time from onset of infectiousness to testing/reporting 
    else
        rep_dist=gampdf(1:25,3.7,3.57);
    end
    elseif scenario==11
    if i<30
        rep_dist=gampdf(1:25,1.85,3.57); % Probability distribution for time from onset of infectiousness to testing/reporting 
    else
        rep_dist=gampdf(1:25,.925,3.57);
    end
    end    
    for j=1:25
        Ntested(i+j,scenario)=max(1,Ntested(i+j,scenario)+round(test_ratio(i,scenario)*NewCase(i,scenario)*rep_dist(j))); % Number tested is equal to the testing ratio on day i times the number of new infections on day i times the reporting delay (pdf that case is reported on day i+j), and at least =1 
        Npositive(i+j,scenario)=min(Ntested(i+j,scenario),Npositive(i+j,scenario)+binornd(NewCase(i,scenario),Ptest_case(i,scenario)*rep_dist(j)));% Number of positive cases is equal to the number of new infections on day i times the probability of testing on day i times the reporting delay (pdf that case is reported on day i+j), and is at most equal to the number tested
    end
end

Ppositive(:,scenario)=Npositive(:,scenario)./Ntested(:,scenario); % Percent of tests positive

% Estimate R0 from the growth rate of the epidemic (days 11-30)
[true_growthrate,~,true_r_stats]=glmfit(1:20,NewCase(21:40,scenario),'Poisson'); % "True" growth rate based on the number of new infections
r_true=true_growthrate(2);

[obs_growthrate,~,obs_r_stats]=glmfit(1:20,Npositive(21:40,scenario),'Poisson'); % "Observed" growth rate based on the number of positive cases
r_obs=obs_growthrate(2);

[test_growthrate,~,test_r_stats]=glmfit(1:20,Ntested(21:40,scenario),'Poisson'); % Growth rate in the number of tests
r_test=test_growthrate(2);

V=1/nu+1/gamma; % Duration of the serial interval (in days)
f=1/(V*nu); % Proportion of the serial interval in the latent period

% Estimate mean value of R0
R0_true(scenario,1)=r_true^2*f*(1-f)*V^2+r_true*V+1;
R0_obs(scenario,1)=r_obs^2*f*(1-f)*V^2+r_obs*V+1;
R0_tests(scenario,1)=r_test^2*f*(1-f)*V^2+r_test*V+1;

% Generate bootstrap samples based on uncertainity in growth rate to
% generate 95% CI for R0 estimates
R0samp_true=zeros(1000,1);
R0samp_obs=zeros(1000,1);
R0samp_tests=zeros(1000,1);
for i=1:1000
    r=normrnd(r_true,true_r_stats.se(2));       
    R0samp_true(i,1)=r^2*f*(1-f)*V^2+r*V+1;
    r=normrnd(r_obs,obs_r_stats.se(2));       
    R0samp_obs(i,1)=r^2*f*(1-f)*V^2+r*V+1;
    r=normrnd(r_test,test_r_stats.se(2));       
    R0samp_tests(i,1)=r^2*f*(1-f)*V^2+r*V+1;
end
R0ci_true(scenario,:)=prctile(R0samp_true,[2.5 97.5]);
R0ci_obs(scenario,:)=prctile(R0samp_obs,[2.5 97.5]);
R0ci_tests(scenario,:)=prctile(R0samp_tests,[2.5 97.5]);

% Estimate Rt using the Wallinga-Teunis method
[Rt_true(:,scenario),Rt_true_adj(:,scenario),Rt_true_ci(:,:,scenario)]=wallinga_teunis(NewCase(1:round(tmax*dt),scenario),[4.79 1.357],100);
[Rt_obs(:,scenario),Rt_obs_adj(:,scenario),Rt_obs_ci(:,:,scenario)]=wallinga_teunis(Npositive(1:round(tmax*dt),scenario),[4.79 1.357],100);
[Rt_obs_full(:,scenario),Rt_obs_full_adj(:,scenario)]=wallinga_teunis(Npositive(:,scenario),[4.79 1.357]);
[Rt_tests(:,scenario),Rt_tests_adj(:,scenario),Rt_tests_ci(:,:,scenario)]=wallinga_teunis(Ntested(1:round(tmax*dt),scenario),[4.79 1.357],100);
[Rt_tests_full(:,scenario),Rt_tests_full_adj(:,scenario)]=wallinga_teunis(Ntested(:,scenario),[4.79 1.357]);

% Estimate the bias in estimates of Rt (as indicated by upper CI<0 or lower CI>0)
bias_upperci(:,scenario)=Rt_obs_ci(:,2,scenario)-Rt_true_adj(:,scenario);
bias_lowerci(:,scenario)=Rt_obs_ci(:,1,scenario)-Rt_true_adj(:,scenario);

end

%% Use EpiEstim to estimate Rt
config = make_config('mean_si',6.5,'std_si',4); %'wind_len',3,
for i=1:11
    Rt_epiestim_true(i)=estimate_R(NewCase(:,i),'parametric_si',config);
    Rt_epiestim_obs(i)=estimate_R(Npositive(1:70,i),'parametric_si',config);
    %Rt_epiestim_obs_full(i)=estimate_R(Npositive(:,i),'parametric_si',config);
end

clear Rt_true_epiestim Rt_true_epiestim_ci Rt_obs_epiestim Rt_obs_epiestim_ci %Rt_obs_full_epiestim Rt_obs_full_epiestim_ci
for i=1:11
    Rt_true_epiestim(:,i)=Rt_epiestim_true(i).median_posterior';
    Rt_true_epiestim_ci(:,:,i)=[Rt_epiestim_true(i).quantile_0025_posterior' Rt_epiestim_true(i).quantile_0975_posterior'];
    Rt_obs_epiestim(:,i)=Rt_epiestim_obs(i).median_posterior';
    Rt_obs_epiestim_ci(:,:,i)=[Rt_epiestim_obs(i).quantile_0025_posterior' Rt_epiestim_obs(i).quantile_0975_posterior'];
    %Rt_obs_full_epiestim(:,i)=Rt_epiestim_obs_full(i).median_posterior';
    %Rt_obs_full_epiestim_ci(:,:,i)=[Rt_epiestim_obs_full(i).quantile_0025_posterior' Rt_epiestim_obs_full(i).quantile_0975_posterior'];
end  

% Estimate the bias in estimates of Rt (as indicated by upper CI<0 or lower CI>0)
for scenario=1:11
    bias_epiestim_upperci(:,scenario)=Rt_obs_epiestim_ci(:,2,scenario)-Rt_true_epiestim(:,scenario);
    bias_epiestim_lowerci(:,scenario)=Rt_obs_epiestim_ci(:,1,scenario)-Rt_true_epiestim(:,scenario);
end

%% PLOT RESULTS %%%

figure
for scenario=1:11
subplot(4,3,scenario)
hold on
bar(NewCase(:,scenario),'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
bar(Ntested(:,scenario),'r','EdgeColor','r') %'LineWidth',1)
bar(Npositive(:,scenario),'b','EdgeColor','b') %'LineWidth',1)
xlim([0 round(tmax*dt)+1])
xlabel('time (days)')
ylabel('# of cases')
title(join({'Scenario',num2str(scenario)}),'FontSize',11)
if scenario==11
    legend('"True" cases','Tests','Observed cases','Location','NW')
end
end

%% Figure 1 main text

figure
for scenario=1:5
subplot(5,4,4*(scenario-1)+1)
%semilogy(cumsum(Npositive(:,scenario)),'b') 
%hold on
%semilogy(cumsum(Ntested(:,scenario)),'r')
%semilogy(cumsum(NewCase(:,scenario)),'--k')
hold on
plot(Ntested(:,scenario),'r','LineWidth',1)
plot(Npositive(:,scenario),'b','LineWidth',1)
if scenario==1
    legend('Tests','Confirmed cases')
end
ylabel('# of cases')
xlim([20 70])
ylim([0 25000])

subplot(5,4,4*(scenario-1)+2)
plot(100*Ppositive(:,scenario),'Color',[.5 0 .5])
ylabel('% positive')
ylim([0 80])
xlim([20 70])
if scenario==1
    legend('Percent positive')
end
if scenario==5
    xlabel('Days since start of epidemic')
end

subplot(5,4,4*(scenario-1)+3)
hold on
%plot(find(Rt_true_adj(:,scenario)),Rt_true_adj(Rt_true_adj(:,scenario)>0,scenario),'k')
%plot(find(Rt_obs_adj(:,scenario)),Rt_obs_adj(Rt_obs_adj(:,scenario)>0,scenario),'Color',[0 .7 0])
%fill([find(Rt_obs_adj(:,scenario)); flipud(find(Rt_obs_adj(:,scenario))); find(Rt_obs_adj(:,scenario),1)],[Rt_obs_ci(find(Rt_obs_adj(:,scenario)),1,scenario); Rt_obs_ci(flipud(find(Rt_obs_adj(:,scenario))),2,scenario); Rt_obs_ci(find(Rt_obs_adj(:,scenario),1),1,scenario)],[.9 1 .9],'EdgeColor',[.9 1 .9])
%fill([find(Rt_true_adj(:,scenario)); flipud(find(Rt_true_adj(:,scenario))); find(Rt_true_adj(:,scenario),1)],[Rt_true_ci(find(Rt_true_adj(:,scenario)),1,scenario); Rt_true_ci(flipud(find(Rt_true_adj(:,scenario))),2,scenario); Rt_true_ci(find(Rt_true_adj(:,scenario),1),1,scenario)],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
%plot(find(Rt_true_adj(:,scenario)),Rt_true_adj(Rt_true_adj(:,scenario)>0,scenario),'k')
%plot(find(Rt_obs_adj(:,scenario)),Rt_obs_adj(Rt_obs_adj(:,scenario)>0,scenario),'Color',[0 .7 0],'LineWidth',1)
%plot(Rt_true(:,scenario),':k')
%plot(Rt_obs(:,scenario),':','Color',[0 .7 0])
plot(Rt_epiestim_true(1).t_end',Rt_true_epiestim(:,scenario),'k')
plot(Rt_epiestim_obs(1).t_end',Rt_obs_epiestim(:,scenario),'Color',[0 .7 0],'LineWidth',1)
fill([Rt_epiestim_obs(1).t_end'; flipud(Rt_epiestim_obs(1).t_end'); Rt_epiestim_obs(1).t_end(1)],[Rt_obs_epiestim_ci(:,1,scenario); flipud(Rt_obs_epiestim_ci(:,2,scenario)); Rt_obs_epiestim_ci(1,1,scenario)],[.9 1 .9],'EdgeColor',[.9 1 .9])
fill([Rt_epiestim_true(1).t_end'; flipud(Rt_epiestim_true(1).t_end'); Rt_epiestim_true(1).t_end(1)],[Rt_true_epiestim_ci(:,1,scenario); flipud(Rt_true_epiestim_ci(:,2,scenario)); Rt_true_epiestim_ci(1,1,scenario)],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
plot(Rt_epiestim_true(1).t_end',Rt_true_epiestim(:,scenario),'k')
plot(Rt_epiestim_obs(1).t_end',Rt_obs_epiestim(:,scenario),'Color',[0 .7 0],'LineWidth',1)
plot((0:round(tmax*dt))',ones(round(tmax*dt)+1,1),'--k')
plot([20 20],[0 6],'k','LineWidth',.75)
ylabel('R_t')
if scenario==1
legend('True','Observed') 
end
xlim([20 67])
ylim([0 6])

subplot(5,4,4*scenario)
hold on
bar(find(bias_lowerci(:,scenario)>0),bias_lowerci(bias_lowerci(:,scenario)>0,scenario),'k')
bar(find(bias_upperci(:,scenario)<0),bias_upperci(bias_upperci(:,scenario)<0,scenario),'k')
xlim([20 67])
ylim([-2 2])
if scenario==1
legend('Bias')
end
end

%%
gtext('A','FontWeight','bold','FontSize',8)
gtext('B','FontWeight','bold','FontSize',8)
gtext('C','FontWeight','bold','FontSize',8)
gtext('D','FontWeight','bold','FontSize',8)
gtext('E','FontWeight','bold','FontSize',8)
%gtext('F','FontWeight','bold')

%% Figure 2 main text

figure
for scenario=6:11
subplot(6,4,4*(scenario-6)+1)
%semilogy(cumsum(Npositive(:,scenario)),'b') 
%hold on
%semilogy(cumsum(Ntested(:,scenario)),'r')
%semilogy(cumsum(NewCase(:,scenario)),'--k')
hold on
plot(Ntested(:,scenario),'r','LineWidth',1)
plot(Npositive(:,scenario),'b','LineWidth',1)
if scenario==6
    legend('Tests','Confirmed cases')
end
ylabel('# of cases')
xlim([20 70])
ylim([0 30000])

subplot(6,4,4*(scenario-6)+2)
plot(100*Ppositive(:,scenario),'Color',[.5 0 .5])
ylabel('% positive')
ylim([0 80])
xlim([20 70])
if scenario==6
    legend('Percent positive')
end
if scenario==11
    xlabel('Days since start of epidemic')
end

subplot(6,4,4*(scenario-6)+3)
hold on
%plot(find(Rt_true_adj(:,scenario)),Rt_true_adj(Rt_true_adj(:,scenario)>0,scenario),'k')
%plot(find(Rt_obs_adj(:,scenario)),Rt_obs_adj(Rt_obs_adj(:,scenario)>0,scenario),'Color',[0 .7 0])
%fill([find(Rt_obs_adj(:,scenario)); flipud(find(Rt_obs_adj(:,scenario))); find(Rt_obs_adj(:,scenario),1)],[Rt_obs_ci(find(Rt_obs_adj(:,scenario)),1,scenario); Rt_obs_ci(flipud(find(Rt_obs_adj(:,scenario))),2,scenario); Rt_obs_ci(find(Rt_obs_adj(:,scenario),1),1,scenario)],[.9 1 .9],'EdgeColor',[.9 1 .9])
%fill([find(Rt_true_adj(:,scenario)); flipud(find(Rt_true_adj(:,scenario))); find(Rt_true_adj(:,scenario),1)],[Rt_true_ci(find(Rt_true_adj(:,scenario)),1,scenario); Rt_true_ci(flipud(find(Rt_true_adj(:,scenario))),2,scenario); Rt_true_ci(find(Rt_true_adj(:,scenario),1),1,scenario)],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
%plot(find(Rt_true_adj(:,scenario)),Rt_true_adj(Rt_true_adj(:,scenario)>0,scenario),'k')
%plot(find(Rt_obs_adj(:,scenario)),Rt_obs_adj(Rt_obs_adj(:,scenario)>0,scenario),'Color',[0 .7 0])
%plot(Rt_true(:,scenario),':k')
%plot(Rt_obs(:,scenario),':','Color',[0 .7 0])
plot(Rt_epiestim_true(1).t_end',Rt_true_epiestim(:,scenario),'k')
plot(Rt_epiestim_obs(1).t_end',Rt_obs_epiestim(:,scenario),'Color',[0 .7 0],'LineWidth',1)
fill([Rt_epiestim_obs(1).t_end'; flipud(Rt_epiestim_obs(1).t_end'); Rt_epiestim_obs(1).t_end(1)],[Rt_obs_epiestim_ci(:,1,scenario); flipud(Rt_obs_epiestim_ci(:,2,scenario)); Rt_obs_epiestim_ci(1,1,scenario)],[.9 1 .9],'EdgeColor',[.9 1 .9])
fill([Rt_epiestim_true(1).t_end'; flipud(Rt_epiestim_true(1).t_end'); Rt_epiestim_true(1).t_end(1)],[Rt_true_epiestim_ci(:,1,scenario); flipud(Rt_true_epiestim_ci(:,2,scenario)); Rt_true_epiestim_ci(1,1,scenario)],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
plot(Rt_epiestim_true(1).t_end',Rt_true_epiestim(:,scenario),'k')
plot(Rt_epiestim_obs(1).t_end',Rt_obs_epiestim(:,scenario),'Color',[0 .7 0],'LineWidth',1)
plot((0:round(tmax*dt))',ones(round(tmax*dt)+1,1),'--k')
plot([20 20],[0 8],'k','LineWidth',.75)
ylabel('R_t')
if scenario==6
legend('True','Observed') 
end
xlim([20 67])
ylim([0 8])

subplot(6,4,4*(scenario-5))
hold on
bar(find(bias_lowerci(:,scenario)>0),bias_lowerci(bias_lowerci(:,scenario)>0,scenario),'k')
bar(find(bias_upperci(:,scenario)<0),bias_upperci(bias_upperci(:,scenario)<0,scenario),'k')
xlim([20 67])
ylim([-5 5])
if scenario==6
legend('Bias')
end

end

%% Estimates of R0

figure
hold on
plot(.9:1:10.9,R0_true,'ok')
plot(1:11,R0_obs,'ob')
plot(1.1:11.1,R0_tests,'or')
plot([.9:1:10.9; .9:1:10.9],R0ci_true','k')
plot([1:11; 1:11],R0ci_obs','b')
plot([1.1:1:11.1; 1.1:1:11.1],R0ci_tests','r')
plot([0 12],[mean(R0_true) mean(R0_true)],'Color',[.7 .7 .7])
legend('True value','Confirmed cases','Number of tests','95% CI')
ylabel('R_0')
xlabel('Scenario')
set(gca,'XTick',1:11)

%% Estimates of Rt based on number of true infections, positive cases, and tests

figure
for scenario=1:11
subplot(4,3,scenario)
hold on
[ax,y1,y2]=plotyy(find(Rt_true_adj(:,scenario)),Rt_true_adj(Rt_true_adj(:,scenario)>0,scenario),1:round(tmax*dt),100*Ppositive(1:round(tmax*dt),scenario));
set(y1(1),'Color','k')
set(y2(1),'Color',[.5 0 .5])
set(ax(1),'YColor','k','YLim',[0 7],'YTick',0:3:6)
set(ax(2),'YColor','k','XLim',[20 60],'YLim',[0 80],'YTick',0:40:80)
plot(find(Rt_obs_adj(:,scenario)),Rt_obs_adj(Rt_obs_adj(:,scenario)>0,scenario),'b')
plot(find(Rt_tests_adj(:,scenario)),Rt_tests_adj(Rt_tests_adj(:,scenario)>0,scenario),'r')
plot(Rt_epiestim_true(1).t_end',Rt_true_epiestim(:,scenario),'--k')
plot(Rt_epiestim_obs(1).t_end',Rt_obs_epiestim(:,scenario),'--b')
if scenario==11
legend('True (WT)','Observed (WT)','Tested (WT)','True (EpiEstim)','Obs (EpiEstim)','% positive')
end
plot((0:round(tmax*dt))',ones(round(tmax*dt)+1,1),'--k')
ylabel('R_t')
ax(2).YLabel.String='% positive';
xlim([20 60])
ylim([0 7])
set(gca,'YTick',0:3:6)
%title(join({'Scenario',num2str(scenario)}),'FontSize',11)
if scenario==11
    xlabel('Days since start of epidemic')
end
end
