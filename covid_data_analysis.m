%% Data cleaning for data downloaded from covidtracking.org 

% Define date range
date=[2020*ones(28,1) 3*ones(28,1) (4:31)'; 2020*ones(30,1) 4*ones(30,1) (1:30)'; 2020*ones(31,1) 5*ones(31,1) (1:31)'; 2020*ones(30,1) 6*ones(30,1) (1:30)'; 2020*ones(31,1) 7*ones(31,1) (1:31)'; 2020*ones(30,1) 8*ones(30,1) (1:30)'];
dateUS=datenum(date);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% US DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Use data import tool to import data from covidtracking.org
npos=flipud(positiveUS);
nneg=flipud(negativeUS);
nhosp=flipud(hospitalizedCumulative);

%% Calculate number of new cases/tests per day and percent positive
ntest=npos+nneg;

newcases=diff(npos); %[npos(1); diff(npos)];
newtests=diff(ntest); %[ntest(1); diff(ntest)];
newhosp=[nhosp(1); max(0,diff(nhosp))];

ppos=newcases./newtests;
phosp=newhosp./newcases;
ppos_mavg=zeros(length(ppos),1);
phosp_mavg=zeros(length(phosp),1);
for i=1:length(newcases)
    ppos_mavg(i,1)=mean(newcases(max(1,i-7):min(length(newcases),i+7)))/mean(newtests(max(1,i-7):min(length(newtests),i+7)));
    phosp_mavg(i,1)=mean(newhosp(max(1,i-7):min(length(newhosp),i+7)))/mean(newcases(max(1,i-7):min(length(newcases),i+7)));
end
%% Calculate the value of R0 for the epidemic using the growth rate ?r? 

%[usgrowthrate,~,us_r_stats]=glmfit((1:21)',npos(1:21,1),'Poisson');
[usgrowthrate,~,us_r_stats]=glmfit((1:21)',newcases(1:21,1),'Poisson');
r_mean=usgrowthrate(2);
%[usgrowthrate_test,~,us_r_test_stats]=glmfit((1:21)',ntest(1:21,1),'Poisson');
[usgrowthrate_test,~,us_r_test_stats]=glmfit((1:21)',newtests(1:21,1),'Poisson');
r_mean_test=usgrowthrate_test(2);

V_mean=6.5; % Duration of the serial interval (in days)
f_mean=2.9/6.5; % Proportion of the serial interval in the latent period

R0mean_US=r_mean^2*f_mean*(1-f_mean)*V_mean^2+r_mean*V_mean+1
R0mean_US_test=r_mean_test^2*f_mean*(1-f_mean)*V_mean^2+r_mean_test*V_mean+1

R0est_US=zeros(1000,1);
for i=1:1000
    r=normrnd(r_mean,us_r_stats.se(2));    
    
    V=V_mean; %gamrnd(4.79,1.357); %
    f=f_mean; %min(.9,gamrnd(3.0585,1.504)/V); %

    R0est_US(i,1)=r^2*f*(1-f)*V^2+r*V+1;
end
R0ci_US=prctile(R0est_US,[2.5 97.5])

R0est_US_test=zeros(1000,1);
for i=1:1000
    r=normrnd(r_mean_test,us_r_test_stats.se(2));    
    
    V=V_mean; %gamrnd(4.79,1.357); %
    f=f_mean; %min(.9,gamrnd(3.0585,1.504)/V); %

    R0est_US_test(i,1)=r^2*f*(1-f)*V^2+r*V+1;
end
R0ci_US_test=prctile(R0est_US_test,[2.5 97.5])

%% Estimate Rt for US data using Wallinga-Teunis
[Rt_US,Rt_US_adj,Rt_US_ci]=wallinga_teunis(newcases,[4.79 1.357],1000);
[Rt_US_tests,Rt_US_tests_adj]=wallinga_teunis(newtests,[4.79 1.357]);

%% Estimate Rt for US data using EpiEstim
config = make_config('mean_si',6.5,'std_si',4); %'wind_len',3,

Rt_epiestim_US=estimate_R(newcases,'parametric_si',config);
Rt_epiestim_UStests=estimate_R(newtests,'parametric_si',config);

Rt_US_epiestim=Rt_epiestim_US.median_posterior';
Rt_US_epiestim_ci=[Rt_epiestim_US.quantile_0025_posterior' Rt_epiestim_US.quantile_0975_posterior'];
Rt_UStests_epiestim=Rt_epiestim_UStests.median_posterior';
Rt_UStests_epiestim_ci=[Rt_epiestim_UStests.quantile_0025_posterior' Rt_epiestim_UStests.quantile_0975_posterior'];

%% Final US figure
figure
subplot(1,4,1)
hold on
[ax,y1,y2]=plotyy(dateUS,newtests,dateUS,newcases);
set(y1,'Color','r','LineWidth',1)
set(y2,'Color','b','LineWidth',1)
set(ax(1),'XLim',[dateUS(1)-3 dateUS(end)+2],'XTick',[dateUS(1)-3; dateUS([59 120],1); dateUS(end)+2],'XTickLabel',{'Mar-01','May-01','Jul-01','Sep-01'},'YColor','k','FontSize',8)
set(ax(2),'XLim',[dateUS(1)-3 dateUS(end)+2],'XTick',[dateUS(1)-3; dateUS([59 120],1); dateUS(end)+2],'XTickLabel',{'Mar-01','May-01','Jul-01','Sep-01'},'YColor','k','FontSize',8)
ax(1).YLabel.String='# of tests';
ax(2).YLabel.String='# of positive cases';
legend('Tests','Positive cases')
    
subplot(1,4,2)
%semilogy(dateUS(1:29),ntest(1:29),'r')
semilogy(dateUS(1:29),newtests(1:29),'r','LineWidth',1)
hold on
%semilogy(dateUS(1:29),npos(1:29),'b')
semilogy(dateUS(1:29),newcases(1:29),'b','LineWidth',1)
semilogy([dateUS(21) dateUS(21)],[1 1e7],'--k')
semilogy(dateUS(1:29),glmval(usgrowthrate,1:29,'log'),'Color',[.7 .7 .7])
semilogy(dateUS(1:29),glmval(usgrowthrate_test,1:29,'log'),'Color',[.7 .7 .7])
datetick('x','mmm-dd-yy')
ylabel('# of tests')
set(gca,'XTick',[dateUS(1)-3; dateUS([12 29],1)],'XTickLabel',{'Mar-01','Mar-15','Apr-01'},'FontSize',8)
xlim([dateUS(1)-3 dateUS(29)+1])
ylim([0 1e6])
legend('Tests','Positive cases','Location','NW')

subplot(1,4,3)
hold on
plot(dateUS,100*ppos,'Color',[.7 .7 .7])
plot(dateUS,100*ppos_mavg,'Color',[.5 0 .5],'LineWidth',1)
ylabel('% positive')
legend('Daily data','Cumulative data')
datetick('x','mmm-dd-yy')
xlim([dateUS(1)-3 dateUS(end)+2])
set(gca,'XTick',[dateUS(1)-3; dateUS([59 120],1); dateUS(end)+2],'XTickLabel',{'Mar-01','May-01','Jul-01','Sep-01'},'FontSize',8)

subplot(1,4,4)
hold on
plot(dateUS(1:end-12,1),Rt_US(1:end-11,1),'Color',[0 .7 0],'LineWidth',1)
plot(dateUS(1:end-1,1),Rt_US(1:end,1),':k')
plot(dateUS(1:end-1,1),Rt_US_adj(1:end,1),'k')
fill([dateUS(1:end-1,1); dateUS(end-1:-1:1,1); dateUS(1)],[Rt_US_ci(1:end,1); Rt_US_ci(end:-1:1,2); Rt_US_ci(1,1)],[.7 1 .7],'EdgeColor',[.7 1 .7])
plot(dateUS(1:end-1,1),Rt_US(1:end,1),':k')
plot(dateUS(1:end-1,1),Rt_US_adj(1:end,1),'k')
plot(dateUS,ones(length(dateUS),1),'--k')
plot(dateUS(1:end-12,1),Rt_US(1:end-11,1),'Color',[0 .7 0],'LineWidth',1)
ylabel('R_t')
datetick('x','mmm-dd-yy')
set(gca,'XTick',[dateUS(1)-3; dateUS([59 120],1)],'XTickLabel',{'Mar-01','May-01','Jul-01'},'FontSize',8)
xlim([dateUS(1)-3 dateUS(end)-3])
ylim([0 5])
legend('Estimated \it{R_t}','Uncorrected','Corrected')
%%
gtext('A)','FontWeight','bold','FontSize',8)
gtext('B)','FontWeight','bold','FontSize',8)
gtext('C)','FontWeight','bold','FontSize',8)
gtext('D)','FontWeight','bold','FontSize',8)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STATE-LEVEL DATA %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Download the historical state data from covidtracking.org API
% Import data using the data import tool 
% Rename first 4 columns "daily_date", "daily_state", "positive", "negative"

x=date(:,1)*10000+date(:,2)*100+date(:,3);
z=categorical(states);
pos_state=zeros(length(dateUS),51);
neg_state=zeros(length(dateUS),51);
for i=1:length(daily_date)
    for j=1:51
        for t=1:length(dateUS)
            if daily_date(i,1)==x(t,1)
                if daily_state(i,1)==z(j)
                    pos_state(t,j)=positive(i,1);
                    neg_state(t,j)=negative(i,1);
                end
            end
        end
    end
end
pos_state(isnan(pos_state))=0; % Replace NaNs with 0
neg_state(isnan(neg_state))=0;

%% Create new summary variables
npos_bystate=pos_state;
nneg_bystate=neg_state;

ntest_bystate=npos_bystate+nneg_bystate;
ppos_cum_bystate=npos_bystate./ntest_bystate;

newcases_bystate=[npos_bystate(1,:); max(0,diff(npos_bystate))];
newtests_bystate=[ntest_bystate(1,:); max(0,diff(ntest_bystate))];

ppos_bystate=zeros(size(newcases_bystate));
for i=1:length(newcases_bystate)
    ppos_bystate(i,:)=mean(newcases_bystate(max(1,i-7):min(length(newcases_bystate),i+7),:))./mean(newtests_bystate(max(1,i-7):min(length(newtests_bystate),i+7),:));
end
%% Calculate the value of R0 for the epidemic using the growth rate ?r? 

V_mean=6.5; % Duration of the serial interval (in days)
f_mean=2.9/6.5; % Proportion of the serial interval in the latent period

for j=1:51
    t0(j)=find(npos_bystate(:,j),1);
    [st_growthrate_cum(:,j),~,st_rcum_stats(j)]=glmfit((t0(j):21)',npos_bystate(t0(j):21,j),'Poisson');
    [st_growthrate(:,j),~,st_r_stats(j)]=glmfit((t0(j):21)',newcases_bystate(t0(j):21,j),'Poisson');
    r_mean(j,1)=st_growthrate(2,j);
    R0mean_st(j,1)=r_mean(j)^2*f_mean*(1-f_mean)*V_mean^2+r_mean(j)*V_mean+1;
    r_cum(j,1)=st_growthrate_cum(2,j);
    R0cum_st(j,1)=r_cum(j)^2*f_mean*(1-f_mean)*V_mean^2+r_cum(j)*V_mean+1;
end

R0est_bystate=zeros(1000,51);
R0cum_bystate=zeros(1000,51);
for i=1:1000
    for j=1:51
        V=V_mean; %gamrnd(4.79,1.357);
        f=f_mean; %min(.9,gamrnd(3.0585,1.504)/V);

        r=normrnd(r_mean(j),st_r_stats(j).se(2)); 
        R0est_bystate(i,j)=r^2*f*(1-f)*V^2+r*V+1;
        
        r=normrnd(r_cum(j),st_rcum_stats(j).se(2)); 
        R0cum_bystate(i,j)=r^2*f*(1-f)*V^2+r*V+1;

    end
end
R0ci_bystate=prctile(R0est_bystate,[2.5 97.5])';
R0ci_cum_bystate=prctile(R0cum_bystate,[2.5 97.5])';

%% Calculate R0 from the growth rate in the number of tests
for j=1:51
    t1(j)=find(ntest_bystate(:,j),1);
    [st_test_growthrate_cum(:,j),~,st_rcum_test_stats(j)]=glmfit((t1(j):21)',ntest_bystate(t1(j):21,j),'Poisson');
    [st_test_growthrate(:,j),~,st_r_test_stats(j)]=glmfit((t1(j):21)',newtests_bystate(t1(j):21,j),'Poisson');
    r_test_mean(j,1)=st_test_growthrate(2,j);
    R0mean_st_test(j,1)=r_test_mean(j)^2*f_mean*(1-f_mean)*V_mean^2+r_test_mean(j)*V_mean+1;
    r_test_cum(j,1)=st_test_growthrate_cum(2,j);
    R0cum_st_test(j,1)=r_test_cum(j)^2*f_mean*(1-f_mean)*V_mean^2+r_test_cum(j)*V_mean+1;
end

R0est_bystate_test=zeros(1000,51);
R0cum_bystate_test=zeros(1000,51);
for i=1:1000
    for j=1:51
        V=V_mean; %gamrnd(4.79,1.357);
        f=f_mean; %min(.9,gamrnd(3.0585,1.504)/V);
        
        r=normrnd(r_test_mean(j),st_r_test_stats(j).se(2));    
        R0est_bystate_test(i,j)=r^2*f*(1-f)*V^2+r*V+1;
        
        r=normrnd(r_test_cum(j),st_rcum_test_stats(j).se(2));    
        R0cum_bystate_test(i,j)=r^2*f*(1-f)*V^2+r*V+1;
    end
end
R0ci_bystate_test=prctile(R0est_bystate_test,[2.5 97.5])';
R0ci_cum_bystate_test=prctile(R0cum_bystate_test,[2.5 97.5])';

%% Calculate Rt for each state using Wallinga-Teunis
Rt_bystate=zeros(size(newcases_bystate,1)-1,51);
Rt_bystate_adj=zeros(size(newcases_bystate,1)-1,51);
Rt_bystate_ci=zeros(size(newcases_bystate,1)-1,2,51);
Rt_bystate_tests=zeros(size(newtests_bystate,1)-1,51);
Rt_bystate_tests_adj=zeros(size(newtests_bystate,1)-1,51);
for j=1:51
    [Rt_bystate(:,j),Rt_bystate_adj(:,j),Rt_bystate_ci(:,:,j)]=wallinga_teunis(newcases_bystate(:,j),[4.79 1.357],1000);
    [Rt_bystate_tests(:,j),Rt_bystate_tests_adj(:,j)]=wallinga_teunis(newtests_bystate(:,j),[4.79 1.357]);
end

%% Estimate Rt for each state using EpiEstim
config = make_config('mean_si',6.5,'std_si',4); %'wind_len',3,

for j=1:51
    Rt_epiestim_bystate(j)=estimate_R(newcases_bystate(:,j),'parametric_si',config);
    Rt_epiestim_bystate_tests(j)=estimate_R(newtests_bystate(:,j),'parametric_si',config);
end
%%
for j=1:51
    Rt_bystate_epiestim(:,j)=Rt_epiestim_bystate(j).median_posterior';
    Rt_bystate_epiestim_ci(:,:,j)=[Rt_epiestim_bystate(j).quantile_0025_posterior' Rt_epiestim_bystate(j).quantile_0975_posterior'];
    Rt_bystate_tests_epiestim(:,j)=Rt_epiestim_bystate_tests(j).median_posterior';
    Rt_bystate_tests_epiestim_ci(:,:,j)=[Rt_epiestim_bystate_tests(j).quantile_0025_posterior' Rt_epiestim_bystate_tests(j).quantile_0975_posterior'];
end

%% Final figure for select states
%Divide into early vs late periods
st=[5 4 19 10 23 35];
n=length(st);
early=29;

for j=1:51
    [Rt_bystate(:,j),Rt_bystate_adj(:,j)]=wallinga_teunis(newcases_bystate(:,j),[4.79 1.357]);
 end

figure
for i=1:n
    subplot(4,n,i)
    %semilogy(dateUS,npos_bystate(:,st(i)),'b')
    %hold on
    %semilogy(dateUS,ntest_bystate(:,st(i)),'r')
    [ax,y1,y2]=plotyy(dateUS(1:early,1),newtests_bystate(1:early,st(i))/1000,dateUS(1:early,1),newcases_bystate(1:early,st(i))/1000);
    set(y1,'Color','r','LineWidth',1)
    set(y2,'Color','b','LineWidth',1)
    set(ax(1),'XLim',[dateUS(1)-3 dateUS(early)+1],'XTick',[dateUS(1)-3; dateUS(29,1)],'XTickLabel',{'Mar-01','Apr-01'},'YColor','k','FontSize',8)
    set(ax(2),'XLim',[dateUS(1)-3 dateUS(early)+1],'XTick',[dateUS(1)-3; dateUS(29,1)],'XTickLabel',{'Mar-01','Apr-01'},'YColor','k','FontSize',8)
    if i==1
        ax(1).YLabel.String='# of tests (x10^3)';
    elseif i==n
        ax(2).YLabel.String={'# positive cases (x10^3)'};
        legend('Tests','Positive cases')
    end
    title(states(st(i)))
      
    subplot(4,n,n+i)
    hold on
    [ax,y1,y2]=plotyy(dateUS(Rt_bystate(1:early,st(i))>0),Rt_bystate(Rt_bystate(1:early,st(i))>0,st(i)),dateUS(1:early,1),100*ppos_bystate(1:early,st(i)));
    set(y1(1),'Color',[0 .7 0],'LineWidth',1)
    set(y2(1),'Color',[.5 0 .5])
    set(ax(1),'XLim',[dateUS(1)-3 dateUS(early)+1],'XTick',[dateUS(1)-3; dateUS(29,1)],'XTickLabel',{'Mar-01','Apr-01'},'YLim',[0 10],'YTick',0:5:10,'YColor','k','FontSize',8)
    set(ax(2),'XLim',[dateUS(1)-3 dateUS(early)+1],'XTick',[dateUS(1)-3; dateUS(29,1)],'XTickLabel',{'Mar-01','Apr-01'},'YLim',[0 100],'YTick',0:50:100,'YColor','k','FontSize',8)
    if i==1
        ax(1).YLabel.String='R_t';
    elseif i==n    
        ax(2).YLabel.String='% positive';
        legend('Estimated R_t','Percent positive')
    end
    fill([dateUS(~isnan(Rt_bystate_ci(1:early,1,st(i))),1); flipud(dateUS(~isnan(Rt_bystate_ci(1:early,1,st(i))),1)); dateUS(find(~isnan(Rt_bystate_ci(1:early,1,st(i))),1))],[Rt_bystate_ci(~isnan(Rt_bystate_ci(1:early,1,st(i))),1,st(i)); flipud(Rt_bystate_ci(~isnan(Rt_bystate_ci(1:early,1,st(i))),2,st(i))); Rt_bystate_ci(find(~isnan(Rt_bystate_ci(1:early,1,st(i))),1),1,st(i))],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
    plot(dateUS(Rt_bystate(1:early,st(i))>0),Rt_bystate(Rt_bystate(1:early,st(i))>0,st(i)),'Color',[0 .7 0],'LineWidth',1)
    plot(dateUS(1:early,1),ones(early,1),'--k')
    
    subplot(4,n,2*n+i)
    %semilogy(dateUS,npos_bystate(:,st(i)),'b')
    %hold on
    %semilogy(dateUS,ntest_bystate(:,st(i)),'r')
    [ax,y1,y2]=plotyy(dateUS(early+1:end,1),newtests_bystate(early+1:end,st(i))/1000,dateUS(early+1:end,1),newcases_bystate(early+1:end,st(i))/1000);
    set(y1,'Color','r','LineWidth',1)
    set(y2,'Color','b','LineWidth',1)
    set(ax(1),'XLim',[dateUS(early) dateUS(end)+1],'XTick',dateUS([59 120],1),'XTickLabel',{'May-01','Jul-01'},'YColor','k','FontSize',8)
    set(ax(2),'XLim',[dateUS(early) dateUS(end)+1],'XTick',dateUS([59 120],1),'XTickLabel',{'May-01','Jul-01'},'YColor','k','FontSize',8)
    if i==1
        ax(1).YLabel.String='# of tests (x10^3)';
    elseif i==n
        ax(2).YLabel.String={'# positive cases (x10^3)'};
    end
    title(states(st(i)))
    
    subplot(4,n,3*n+i)
    hold on
    [ax,y1,y2]=plotyy(dateUS(Rt_bystate(1:end-12,st(i))>0),Rt_bystate(Rt_bystate(1:end-12,st(i))>0,st(i)),dateUS(early+1:end,1),100*ppos_bystate(early+1:end,st(i)));
    set(y1(1),'Color',[0 .7 0],'LineWidth',1)
    set(y2(1),'Color',[.5 0 .5])
    set(ax(1),'XLim',[dateUS(early) dateUS(end)+1],'XTick',dateUS([59 120],1),'XTickLabel',{'May-01','Jul-01'},'YLim',[0 2],'YTick',0:2,'YColor','k','FontSize',8)
    set(ax(2),'XLim',[dateUS(early) dateUS(end)+1],'XTick',dateUS([59 120],1),'XTickLabel',{'May-01','Jul-01'},'YLim',[0 40],'YTick',0:20:40,'YColor','k','FontSize',8)
    if i==1
        ax(1).YLabel.String='R_t';
    elseif i==n    
        ax(2).YLabel.String='% positive';
    end
    fill([dateUS(~isnan(Rt_bystate_ci(1:end,1,st(i))),1); flipud(dateUS(~isnan(Rt_bystate_ci(1:end,1,st(i))),1)); dateUS(find(~isnan(Rt_bystate_ci(1:end,1,st(i))),1))],[Rt_bystate_ci(~isnan(Rt_bystate_ci(1:end,1,st(i))),1,st(i)); flipud(Rt_bystate_ci(~isnan(Rt_bystate_ci(1:end,1,st(i))),2,st(i))); Rt_bystate_ci(find(~isnan(Rt_bystate_ci(1:end,1,st(i))),1),1,st(i))],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
    plot(dateUS(Rt_bystate_adj(1:end,st(i))>0),Rt_bystate_adj(Rt_bystate_adj(1:end,st(i))>0,st(i)),'k')
    plot(dateUS(end-13:end-1),Rt_bystate(end-12:end,st(i)),':k')
    plot(dateUS(Rt_bystate(1:end-12,st(i))>0),Rt_bystate(Rt_bystate(1:end-12,st(i))>0,st(i)),'Color',[0 .7 0],'LineWidth',1)
    plot(dateUS(early+1:end),ones(length(dateUS)-early,1),'--k')
    plot([dateUS(early) dateUS(early)],[0 2],'k')
    
end

%%
gtext('A)','FontWeight','bold','FontSize',8)
gtext('B)','FontWeight','bold','FontSize',8)
gtext('C)','FontWeight','bold','FontSize',8)
gtext('D)','FontWeight','bold','FontSize',8)

%% Per capita testing rates and growth rates for select states
for i=1:51
    [st_test_r(:,i),~,st_test_rstats(i)]=glmfit((1:91)',10000*newtests_bystate(29:119,i)./pop_bystate(i),'Poisson');
    [st_test_rinit(:,i),~,st_test_rinitstats(i)]=glmfit((t1(i):21)',10000*newtests_bystate(t1(i):21,i)./pop_bystate(i),'Poisson');
end

figure
for i=1:6
subplot(2,6,i)
hold on
plot(dateUS(1:21),glmval(st_test_rinit,(1:21)','log'),'Color',[.7 .7 .7])
plot(dateUS(1:28),10000*newtests_bystate(1:28,st(i))./pop_bystate(st(i)),'r')
plot(dateUS(1:21),glmval(st_test_rinit(:,st(i)),(1:21)','log'),'k')
datetick('x','mmm-dd')
xlim([dateUS(1)-3 dateUS(29)])
ylim([0 15])
if i==1
ylabel('# of tests (per 10K)')
end
title(states(st(i)))

subplot(2,6,i+6)
hold on
plot(dateUS(29:119),glmval(st_test_r,(1:91)','log'),'Color',[.7 .7 .7])
plot(dateUS(29:end),10000*newtests_bystate(29:end,st(i))./pop_bystate(st(i)),'r')
plot(dateUS(29:119),glmval(st_test_r(:,st(i)),(1:91)','log'),'k')
datetick('x','mmm-dd')
xlim([dateUS(28) dateUS(end)+1])
set(gca,'XTick',dateUS([29 90 151]),'XTickLabel',{'Apr-01','Jun-01','Aug-01'})
ylim([0 50])
if i==1
ylabel('# of tests (per 10K)')
end
title(states(st(i)))
end
%%
gtext('A)','FontWeight','bold')
gtext('B)','FontWeight','bold')

%% Plot Rt for all 51 US states and territories
for j=1:51
    [Rt_bystate(:,j),Rt_bystate_adj(:,j)]=wallinga_teunis(newcases_bystate(:,j),[4.79 1.357]);
    [Rt_bystate_tests(:,j),Rt_bystate_tests_adj(:,j)]=wallinga_teunis(newtests_bystate(:,j),[4.79 1.357]);
end
%%
figure
for i=1:51
    subplot(9,6,i)
    hold on
    %[ax,y1,y2]=plotyy(dateUS(Rt_bystate(1:end-12,i)>0),Rt_bystate(Rt_bystate(1:end-12,i)>0,i),dateUS,100*ppos_bystate(:,i));
    [ax,y1,y2]=plotyy(dateUS(Rt_epiestim_bystate(i).t_end'),Rt_bystate_epiestim(:,i),dateUS,100*ppos_bystate(:,i));
    set(y1(1),'Color',[0 .7 0],'LineWidth',1)
    set(y2(1),'Color',[.5 0 .5])
    set(ax(1),'YColor','k','FontSize',8)
    set(ax(2),'YColor','k','XLim',[dateUS(1)-3 dateUS(end)+1],'YLim',[0 100],'YTick',0:50:100,'FontSize',8)
    %plot(dateUS(Rt_bystate_tests_adj(1:end-1,i)>0),Rt_bystate_tests_adj(Rt_bystate_tests_adj(1:end-1,i)>0,i),'r')
    if i==51
        legend('R_t','Percent positive')
    end
    %fill([dateUS(~isnan(Rt_bystate_ci(:,1,i)),1); flipud(dateUS(~isnan(Rt_bystate_ci(:,1,i)),1)); dateUS(find(~isnan(Rt_bystate_ci(:,1,i)),1))],[Rt_bystate_ci(~isnan(Rt_bystate_ci(:,1,i)),1,i); flipud(Rt_bystate_ci(~isnan(Rt_bystate_ci(:,1,i)),2,i)); Rt_bystate_ci(find(~isnan(Rt_bystate_ci(:,1,i)),1),1,i)],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
    %plot(dateUS(Rt_bystate_tests_adj(1:end-1,i)>0),Rt_bystate_tests_adj(Rt_bystate_tests_adj(1:end-1,i)>0,i),'r')
    %plot(dateUS(Rt_bystate_adj(1:end-1,i)>0),Rt_bystate_adj(Rt_bystate_adj(1:end-1,i)>0,i),'k')
    %plot(dateUS(Rt_bystate(1:end-12,i)>0),Rt_bystate(Rt_bystate(1:end-12,i)>0,i),'Color',[0 .7 0],'LineWidth',1)
    fill([dateUS(Rt_epiestim_bystate(i).t_end'); dateUS(flipud(Rt_epiestim_bystate(i).t_end')); dateUS(Rt_epiestim_bystate(i).t_end(1))],[Rt_bystate_epiestim_ci(:,1,i); flipud(Rt_bystate_epiestim_ci(:,2,i)); Rt_bystate_epiestim_ci(1,1,i)],[.9 1 .9],'EdgeColor',[.9 1 .9])
    plot(dateUS(Rt_epiestim_bystate(i).t_end'),Rt_bystate_epiestim(:,i),'Color',[0 .7 0],'Linewidth',1)
    plot(dateUS,ones(length(dateUS),1),'--k')
    if rem(i,6)==1
    ylabel('R_t')
    end
    if rem(i,6)==0
    ax(2).YLabel.String='% positive';
    end
    datetick('x','mmm-dd')
    xlim([dateUS(1)-3 dateUS(end)+2])
    ylim([0 15])
    set(gca,'YTick',0:5:15,'XTick',[dateUS(1)-3 dateUS(90) dateUS(end)+2],'XTickLabel',{'Mar1','Jun1','Sep1'},'FontSize',8)
    title(states(i)) 
end
    

subplot(9,6,54)
hold on
[ax,y1,y2]=plotyy(dateUS(1:end-13,1),Rt_US(1:end-12,1),dateUS,ppos_mavg);
set(y1(1),'Color',[0 .7 0],'LineWidth',1)
set(y2(1),'Color',[.5 0 .5])
set(ax(1),'YColor','k','FontSize',8)
set(ax(2),'YColor','k','XLim',[dateUS(1)-3 dateUS(end)+1],'YLim',[0 100],'YTick',0:50:100,'FontSize',8)
fill([dateUS(1:end-1,1); dateUS(end-1:-1:1,1); dateUS(1)],[Rt_US_ci(1:end,1); Rt_US_ci(end:-1:1,2); Rt_US_ci(1,1)],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
%plot(dateUS(1:end-1,1),Rt_US_tests_adj,'r')
plot(dateUS(1:end-1,1),Rt_US,':k')
plot(dateUS(1:end-1,1),Rt_US_adj,'k')
plot(dateUS(1:end-13,1),Rt_US(1:end-12,1),'Color',[0 .7 0],'LineWidth',1)
plot(dateUS,ones(length(dateUS),1),'--k')
ylabel('R_t')
ax(2).YLabel.String='% positive';
datetick('x','mmm')
xlim([dateUS(1)-3 dateUS(end)+2])
ylim([0 15])
set(gca,'YTick',0:5:15,'XTick',[dateUS(1)-3 dateUS(90) dateUS(end)+2],'XTickLabel',{'Mar1','Jun1','Sep1'},'FontSize',8)
title('US')
