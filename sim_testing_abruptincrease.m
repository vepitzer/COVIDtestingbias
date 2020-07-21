%%% PARAMETER VALUES %%%

tmax=160; %maximum number of time-steps for the simulation
dt=0.5; %time step = 1/20 day (for exponential distributed rates, approximate method)

beta1=[.85*ones(round(tmax*dt)-20,1); .5*ones(20,1)];  %transmission parameter
nu=1/2.9;  %rate of leaving latent period 
gamma=1/3.7;  %rate of recovery

%%% DEFINE TESTING SCENARIOS

rep_dist_pars=gampdf(1:25,1.85,3.57); % Default probability distribution for time from onset of infectiousness to testing/reporting 

Ptest_case=.1*ones(round(tmax*dt),11); % Probability that a "true" case is tested (default scenario)
Ptest_case(:,1)=[.05*ones(round(tmax*dt)-50,1); .15*ones(50,1)];
Ptest_case(:,2)=[.05*ones(round(tmax*dt)-50,1); .25*ones(50,1)];
Ptest_case(:,3)=[.05*ones(round(tmax*dt)-50,1); .35*ones(50,1)];
Ptest_case(:,4)=[.05*ones(round(tmax*dt)-50,1); .45*ones(50,1)];
Ptest_case(:,5)=[.05*ones(round(tmax*dt)-50,1); .55*ones(50,1)];

Ptest_case(:,11)=[.05*ones(round(tmax*dt)-50,1); .15*ones(50,1)];
Ptest_case(:,12)=[.05*ones(round(tmax*dt)-50,1); .25*ones(50,1)];
Ptest_case(:,13)=[.05*ones(round(tmax*dt)-50,1); .35*ones(50,1)];
Ptest_case(:,14)=[.05*ones(round(tmax*dt)-50,1); .45*ones(50,1)];
Ptest_case(:,15)=[.05*ones(round(tmax*dt)-50,1); .55*ones(50,1)];

test_ratio=ones(round(tmax*dt),11); % Ratio of the number tested to the number of "true" cases (default scenario)
test_ratio(:,6)=[.2*ones(round(tmax*dt)-50,1); .4*ones(50,1)]; % Ratio of the number tested to the number of "true" cases
test_ratio(:,7)=[.2*ones(round(tmax*dt)-50,1); .6*ones(50,1)]; % Ratio of the number tested to the number of "true" cases
test_ratio(:,8)=[.2*ones(round(tmax*dt)-50,1); .8*ones(50,1)]; % Ratio of the number tested to the number of "true" cases
test_ratio(:,9)=[.2*ones(round(tmax*dt)-50,1); 1*ones(50,1)]; % Ratio of the number tested to the number of "true" cases
test_ratio(:,10)=[.2*ones(round(tmax*dt)-50,1); 2*ones(50,1)]; % Ratio of the number tested to the number of "true" cases

test_ratio(:,11)=[.2*ones(round(tmax*dt)-50,1); .4*ones(50,1)]; % Ratio of the number tested to the number of "true" cases
test_ratio(:,12)=[.2*ones(round(tmax*dt)-50,1); .6*ones(50,1)]; % Ratio of the number tested to the number of "true" cases
test_ratio(:,13)=[.2*ones(round(tmax*dt)-50,1); .8*ones(50,1)]; % Ratio of the number tested to the number of "true" cases
test_ratio(:,14)=[.2*ones(round(tmax*dt)-50,1); 1*ones(50,1)]; % Ratio of the number tested to the number of "true" cases
test_ratio(:,15)=[.2*ones(round(tmax*dt)-50,1); 2*ones(50,1)]; % Ratio of the number tested to the number of "true" cases


%%% INITIALIZE OUTPUT VARIABLES %%%
NewCase_delta=zeros(round(tmax*dt),15);
Ntested_delta=ones(round(tmax*dt)+25,15);
Npositive_delta=zeros(round(tmax*dt)+25,15);
Ppositive_delta=zeros(round(tmax*dt)+25,15);

R0_true_delta=zeros(15,1);
R0_obs_delta=zeros(15,1);
R0_tests_delta=zeros(15,1);
R0ci_true_delta=zeros(15,2);
R0ci_obs_delta=zeros(15,2);
R0ci_tests_delta=zeros(15,2);

Rt_true_delta=zeros(round(tmax*dt)-1,15);
Rt_true_delta_ci=zeros(round(tmax*dt)-1,2,15);
Rt_obs_delta=zeros(round(tmax*dt)-1,15);
Rt_obs_delta_ci=zeros(round(tmax*dt)-1,2,15);
Rt_tests_delta=zeros(round(tmax*dt)-1,15);
Rt_tests_delta_ci=zeros(round(tmax*dt)-1,2,15);
Rt_bias_delta=zeros(round(tmax*dt)-1,16);

for scenario=1:15

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

tstep=round(1/dt); %Number of time steps per day (= 1/dt)
NewCase_delta(1,scenario)=10; 
for i=2:round(tmax*dt)
    NewCase_delta(i,scenario)=sum(Infectious(tstep*(i-1)+1:tstep*i,1)); %Number of new cases per day is the sum of individuals with symptom onset occurring during the 20 time steps that make up that day
end


%%% Incorporate testing practices and delays

for i=1:round(tmax*dt)  
    Ntested_delta(1,scenario)=1;
    for j=1:25
        Ntested_delta(i+j,scenario)=max(1,Ntested_delta(i+j,scenario)+round(test_ratio(i,scenario)*NewCase_delta(i,scenario)*rep_dist(j)));
        Npositive_delta(i+j,scenario)=min(Ntested_delta(i+j,scenario),Npositive_delta(i+j,scenario)+binornd(NewCase_delta(i,scenario),Ptest_case(i,scenario)*rep_dist(j)));
    end
end

Ppositive_delta(:,scenario)=Npositive_delta(:,scenario)./Ntested_delta(:,scenario);

% Estimate R0 from the growth rate of the epidemic (days 11-30)
[true_growthrate,~,true_r_stats]=glmfit(1:30,NewCase_delta(1:30,scenario),'Poisson');
r_true=true_growthrate(2);
[obs_growthrate,~,obs_r_stats]=glmfit(1:20,Npositive_delta(21:40,scenario),'Poisson');
r_obs=obs_growthrate(2);
[test_growthrate,~,test_r_stats]=glmfit(1:20,Ntested_delta(21:40,scenario),'Poisson');
r_test=test_growthrate(2);

V=1/nu+1/gamma; % Duration of the serial interval (in days)
f=1/(V*nu); % Proportion of the serial interval in the latent period

R0_true_delta(scenario,1)=r_true^2*f*(1-f)*V^2+r_true*V+1;
R0_obs_delta(scenario,1)=r_obs^2*f*(1-f)*V^2+r_obs*V+1;
R0_tests_delta(scenario,1)=r_test^2*f*(1-f)*V^2+r_test*V+1;

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
R0ci_true_delta(scenario,:)=prctile(R0samp_true,[2.5 97.5]);
R0ci_obs_delta(scenario,:)=prctile(R0samp_obs,[2.5 97.5]);
R0ci_tests_delta(scenario,:)=prctile(R0samp_tests,[2.5 97.5]);


% Estimate Rt using the Wallinga-Teunis method
[~,Rt_true_delta(:,scenario),Rt_true_delta_ci(:,:,scenario)]=wallinga_teunis(NewCase_delta(1:round(tmax*dt),scenario),[4.79 1.357],100);
[~,Rt_obs_delta(:,scenario),Rt_obs_delta_ci(:,:,scenario)]=wallinga_teunis(Npositive_delta(1:round(tmax*dt),scenario),[4.79 1.357],100);
[~,Rt_tests_delta(:,scenario),Rt_tests_delta_ci(:,:,scenario)]=wallinga_teunis(Ntested_delta(1:round(tmax*dt),scenario),[4.79 1.357],100);

Rt_bias_delta(:,scenario)=max(0,Rt_obs_delta_ci(:,1,scenario)-Rt_true_delta(:,scenario));
end


%% PLOT RESULTS %%%

col=colormap(lines(5));

figure
for scenario=1:3
subplot(3,4,4*(scenario-1)+1)
hold on
for i=1:5
plot(Ntested_delta(:,(scenario-1)*5+i),'--','Color',col(i,:),'LineWidth',1)
plot(Npositive_delta(:,(scenario-1)*5+i),'Color',col(i,:),'LineWidth',1)
end
ylabel('# of cases')
%ylim([0 1e6])
xlim([20 70])

subplot(3,4,4*(scenario-1)+2)
plot(100*Ppositive_delta(:,(scenario-1)*5+1:scenario*5))
ylabel('% positive')
ylim([0 100])
xlim([20 70])
if scenario==3
    xlabel('Days since start of epidemic')
end

subplot(3,4,4*(scenario-1)+3)
hold on
plot(Rt_true_delta(:,(scenario-1)*5+1:scenario*5),'k')
plot(Rt_obs_delta(:,(scenario-1)*5+1:scenario*5))
plot((0:round(tmax*dt))',ones(round(tmax*dt)+1,1),'--k')
plot([20 20],[0 6],'k','LineWidth',.75)
ylabel('R_t')
xlim([20 70])
ylim([0 15])

subplot(3,4,4*scenario)
hold on
plot(Rt_bias_delta(:,(scenario-1)*5+1:scenario*5))
xlim([20 70])
ylim([-1 8])
if scenario==1
    legend('10% increase','20% increase','30% increase','40% increase','50% increase')
elseif scenario==2
    legend('2x increase','3x increase','4x increase','5x increase','10x increase')
elseif scenario==3
    legend('10%/2x increase','20%/3x increase','30%/4x increase','40%/5x increase','50%/10x increase')
end
end



%%
%gtext('A','FontWeight','bold')
%gtext('B','FontWeight','bold')
%gtext('C','FontWeight','bold')
%gtext('D','FontWeight','bold')
%gtext('E','FontWeight','bold')
%gtext('F','FontWeight','bold')


%%
figure
hold on
plot(.9:1:14.9,R0_true_delta,'ok')
plot(1:15,R0_obs_delta,'ob')
plot(1.1:15.1,R0_tests_delta,'or')
plot([.9:1:14.9; .9:1:14.9],R0ci_true_delta','k')
plot([1:15; 1:15],R0ci_obs_delta','b')
plot([1.1:1:15.1; 1.1:1:15.1],R0ci_tests_delta','r')
plot([0 16],[mean(R0_true_delta) mean(R0_true_delta)],'Color',[.7 .7 .7])
legend('True value','Confirmed cases','Number of tests','95% CI')
ylabel('R_0')
xlabel('Scenario')
set(gca,'XTick',1:15)
