%%% PARAMETER VALUES %%%

tmax=1400; %maximum number of time-steps for the simulation
dt=0.05; %time step = 1/20 day (for exponential distributed rates, approximate method)

beta1=[.85*ones(round(tmax*dt)-20,1); .3*ones(20,1)];  %transmission parameter
nu=1/3.7;  %rate of leaving latent period 
gamma=1/3.5;  %rate of recovery

%%% DEFINE TESTING SCENARIOS

rep_dist_pars=gampdf(1:25,1.85,3.57); % Default probability distribution for time from onset of infectiousness to testing/reporting 

Ptest_case=zeros(round(tmax*dt),11);
Ptest_case(:,1)=.1*ones(round(tmax*dt),1); % Probability that a "true" case is tested
Ptest_case(:,2)=[.05*ones(round(tmax*dt)-51,1); (.05:.0025:.15)'; .15*ones(10,1)];
Ptest_case(:,3)=[.15*ones(round(tmax*dt)-51,1); (.15:-.0025:.05)'; .05*ones(10,1)];
Ptest_case(:,4:5)=.1*ones(round(tmax*dt),2); % Probability that a "true" case is tested
Ptest_case(:,6)=[.05*ones(round(tmax*dt)-30,1); .15*ones(30,1)];
Ptest_case(:,7)=[.15*ones(round(tmax*dt)-30,1); .05*ones(30,1)];
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
Ntested=zeros(round(tmax*dt)+25,11);
Npositive=zeros(round(tmax*dt)+25,11);

R0_true=zeros(11,1);
R0_obs=zeros(11,1);
R0_tests=zeros(11,1);
R0ci_true=zeros(11,2);
R0ci_obs=zeros(11,2);
R0ci_tests=zeros(11,2);
Rt_true=zeros(round(tmax*dt),11);
Rt_true_adj=zeros(round(tmax*dt),11);
Rt_obs=zeros(round(tmax*dt),11);
Rt_obs_adj=zeros(round(tmax*dt),11);
Rt_obs_full=zeros(round(tmax*dt)+25,11);
Rt_obs_full_adj=zeros(round(tmax*dt)+25,11);
Rt_tests=zeros(round(tmax*dt),11);
Rt_tests_adj=zeros(round(tmax*dt),11);
Rt_tests_full=zeros(round(tmax*dt)+25,11);
Rt_tests_full_adj=zeros(round(tmax*dt)+25,11);

for scenario=1:11

%%% INITIALIZE STATE VARIABLES %%%
S=zeros(tmax,1);
E=zeros(tmax,1);
I=zeros(tmax,1);
R=zeros(tmax,1);
N=zeros(tmax,1);

Infect=zeros(tmax,1);
Symp=zeros(tmax,1);
Recover=zeros(tmax,1);

N(1)=1000000;  %initial number of people 

E(1)=0;
I(1)=10;
R(1)=0;
S(1)=N(1)-E(1)-I(1)-R(1); 
Symp(1)=I(1);

%%% START OF LOOP %%%
time=zeros(tmax,1);
for t=2:tmax
    
%%% RATES %%%
%%% Assuming exponentially distributed latent/infectious period
%%% (Approximate method)
Infect(t)=poissrnd(S(t-1)*beta1(ceil(t*dt))*I(t-1)/N(t-1)*dt); %Number of events is Poisson distributed if we assume they occur at a constant rate over short time step
Symp(t)=poissrnd(nu*E(t-1)*dt);
Recover(t)=poissrnd(gamma*I(t-1)*dt);


%%% TRANSMISSION DYNAMICS %%%

%People
S(t)=max(S(t-1)-Infect(t),0);
E(t)=max(E(t-1)+Infect(t)-Symp(t),0);
I(t)=max(I(t-1)+Symp(t)-Recover(t),0);
R(t)=max(R(t-1)+Recover(t),0);
N(t)=S(t)+E(t)+I(t)+R(t);

%%% END OF LOOP
%time(t)=time(t-1)+tte;
if I(t)==0, break, end

end
tf=t;

tstep=round(1/dt); %Number of time steps per day (= 1/.05 = 20)
NewCase(1,scenario)=10; %1; %
for i=2:round(tmax*dt)
    NewCase(i,scenario)=sum(Symp(tstep*(i-1)+1:tstep*i,1)); %Number of new cases per day is the sum of individuals with symptom onset occurring during the 20 time steps that make up that day
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
        Ntested(i+j,scenario)=Ntested(i+j,scenario)+round(test_ratio(i,scenario)*NewCase(i,scenario)*rep_dist(j));
        Npositive(i+j,scenario)=Npositive(i+j,scenario)+binornd(NewCase(i,scenario),Ptest_case(i,scenario)*rep_dist(j));
    end
end

% Estimate R0 from the growth rate of the epidemic (days 11-30)
[true_growthrate,~,true_r_stats]=glmfit(1:20,cumsum(NewCase(21:40,scenario)),'Poisson');
r_true=true_growthrate(2);
[obs_growthrate,~,obs_r_stats]=glmfit(1:20,cumsum(Npositive(21:40,scenario)),'Poisson');
r_obs=obs_growthrate(2);
[test_growthrate,~,test_r_stats]=glmfit(1:20,cumsum(Ntested(21:40,scenario)),'Poisson');
r_test=test_growthrate(2);

V=1/nu+1/gamma; % Duration of the serial interval (in days)
f=1/(V*nu); % Proportion of the serial interval in the latent period

R0_true(scenario,1)=r_true^2*f*(1-f)*V^2+r_true*V+1;
R0_obs(scenario,1)=r_obs^2*f*(1-f)*V^2+r_obs*V+1;
R0_tests(scenario,1)=r_test^2*f*(1-f)*V^2+r_test*V+1;

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
[Rt_true(:,scenario),Rt_true_adj(:,scenario)]=wallinga_teunis(NewCase(:,scenario),[4.79 1.357]);
[Rt_obs(:,scenario),Rt_obs_adj(:,scenario)]=wallinga_teunis(Npositive(1:round(tmax*dt),scenario),[4.79 1.357]);
[Rt_obs_full(:,scenario),Rt_obs_full_adj(:,scenario)]=wallinga_teunis(Npositive(:,scenario),[4.79 1.357]);
[Rt_tests(:,scenario),Rt_tests_adj(:,scenario)]=wallinga_teunis(Ntested(1:round(tmax*dt),scenario),[4.79 1.357]);
[Rt_tests_full(:,scenario),Rt_tests_full_adj(:,scenario)]=wallinga_teunis(Ntested(:,scenario),[4.79 1.357]);

end

Ppositive=Npositive./Ntested;


%% PLOT RESULTS %%%
%figure
%hold on
%bar(NewCase(:,1))
%bar(Ntested(:,1),'r')
%bar(Npositive(:,1),'b')
%legend('"True" cases','Tests','Observed cases','Location','NW')
%xlim([0 round(tmax*dt)+30])
%xlabel('time (days)')
%ylabel('# of cases')

%% PLOT RESULTS %%%
figure
for scenario=1:11
subplot(4,3,scenario)
hold on
bar(NewCase(:,scenario),'FaceColor',[.7 .7 .7],'EdgeColor',[.7 .7 .7])
bar(Ntested(:,scenario),'r','EdgeColor','r')
bar(Npositive(:,scenario),'b','EdgeColor','b')
xlim([0 round(tmax*dt)+1])
xlabel('time (days)')
ylabel('# of cases')
title(join({'Scenario',num2str(scenario)}),'FontSize',11)
if scenario==11
    legend('"True" cases','Tests','Observed cases','Location','NW')
end
end

%%
figure
for scenario=1:5
subplot(5,3,3*(scenario-1)+1)
semilogy(cumsum(Npositive(:,scenario)),'b') 
hold on
semilogy(cumsum(Ntested(:,scenario)),'r')
semilogy(cumsum(NewCase(:,scenario)),'--k')
if scenario==1
legend('Confirmed cases','Tests','"True" cases','Location','NW')
end
ylabel('# of cases')
ylim([0 1e6])
xlim([20 60])

subplot(5,3,3*(scenario-1)+2)
plot(100*Ppositive(:,scenario),'Color',[.5 0 .5])
ylabel('% positive')
ylim([0 80])
xlim([20 60])
if scenario==1
    legend('Percent positive')
end
if scenario==5
    xlabel('Days since start of epidemic')
end

subplot(5,3,3*scenario)
hold on
plot(find(Rt_true_adj(:,scenario)),Rt_true_adj(Rt_true_adj(:,scenario)>0,scenario),'k')
plot(find(Rt_obs_adj(:,scenario)),Rt_obs_adj(Rt_obs_adj(:,scenario)>0,scenario),'Color',[0 .7 0])
%plot(find(Rt_tests_adj(:,scenario)),Rt_tests_adj(Rt_tests_adj(:,scenario)>0,scenario),'r')
plot(Rt_true(:,scenario),':k')
plot(Rt_obs(:,scenario),':','Color',[0 .7 0])
%plot(Rt_tests(:,scenario),':r')
plot((0:round(tmax*dt))',ones(round(tmax*dt)+1,1),'--k')
ylabel('R_t')
if scenario==1
legend('True','Observed') %,'Tested')
end
xlim([20 60])
ylim([0 6])
end

%%
gtext('A','FontWeight','bold')
gtext('B','FontWeight','bold')
gtext('C','FontWeight','bold')
gtext('D','FontWeight','bold')
gtext('E','FontWeight','bold')
%gtext('F','FontWeight','bold')

%%
figure
for scenario=6:11
subplot(6,3,3*(scenario-6)+1)
semilogy(cumsum(Npositive(:,scenario)),'b') 
hold on
semilogy(cumsum(Ntested(:,scenario)),'r')
semilogy(cumsum(NewCase(:,scenario)),'--k')
if scenario==6
legend('Confirmed cases','Tests','"True" cases','Location','NW')
end
ylabel('# of cases')
ylim([0 1e6])
xlim([20 60])

subplot(6,3,3*(scenario-6)+2)
plot(100*Ppositive(:,scenario),'Color',[.5 0 .5])
ylabel('% positive')
ylim([0 80])
xlim([20 60])
if scenario==6
    legend('Percent positive')
end
if scenario==11
    xlabel('Days since start of epidemic')
end

subplot(6,3,3*(scenario-5))
hold on
plot(find(Rt_true_adj(:,scenario)),Rt_true_adj(Rt_true_adj(:,scenario)>0,scenario),'k')
plot(find(Rt_obs_adj(:,scenario)),Rt_obs_adj(Rt_obs_adj(:,scenario)>0,scenario),'Color',[0 .7 0])
%plot(find(Rt_tests_adj(:,scenario)),Rt_tests_adj(Rt_tests_adj(:,scenario)>0,scenario),'r')
plot(Rt_true(:,scenario),':k')
plot(Rt_obs(:,scenario),':','Color',[0 .7 0])
%plot(Rt_tests(:,scenario),':r')
plot((0:round(tmax*dt))',ones(round(tmax*dt)+1,1),'--k')
ylabel('R_t')
if scenario==6
legend('True','Observed') %,'Tested')
end
xlim([20 60])
ylim([0 6])
end

%%
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

%%
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
if scenario==11
legend('True','Observed','Tested','% positive')
end
plot((0:round(tmax*dt))',ones(round(tmax*dt)+1,1),'--k')
ylabel('R_t')
ax(2).YLabel.String='% positive';
xlim([20 60])
ylim([0 7])
set(gca,'YTick',0:3:6)
title(join({'Scenario',num2str(scenario)}),'FontSize',11)
if scenario==11
    xlabel('Days since start of epidemic')
end
end