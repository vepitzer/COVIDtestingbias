%% Figure 1 main text

figure
for scenario=1:5
subplot(5,4,4*(scenario-1)+1)
%semilogy(cumsum(Npositive(:,scenario)),'b') 
%hold on
%semilogy(cumsum(Ntested(:,scenario)),'r')
%semilogy(cumsum(NewCase(:,scenario)),'--k')
hold on
plot(Ntested(:,scenario),'--k','LineWidth',1)
plot(Npositive(:,scenario),'k','LineWidth',1)
if scenario==1
    legend('Tests','Confirmed cases')
end
ylabel('# of cases')
xlim([20 70])
ylim([0 25000])

subplot(5,4,4*(scenario-1)+2)
plot(100*Ppositive(:,scenario),'k')
ylabel('% positive')
ylim([0 80])
xlim([20 70])
if scenario==1
    legend('%_{pos}')
end
if scenario==5
    xlabel('Days since start of epidemic')
end

subplot(5,4,4*(scenario-1)+3)
hold on
plot(Rt_epiestim_true(1).t_end',Rt_true_epiestim(:,scenario),'k','LineWidth',1)
plot(Rt_epiestim_obs(1).t_end',Rt_obs_epiestim(:,scenario),'--k','LineWidth',1)
fill([Rt_epiestim_obs(1).t_end'; flipud(Rt_epiestim_obs(1).t_end'); Rt_epiestim_obs(1).t_end(1)],[Rt_obs_epiestim_ci(:,1,scenario); flipud(Rt_obs_epiestim_ci(:,2,scenario)); Rt_obs_epiestim_ci(1,1,scenario)],[.7 .7 .7],'EdgeColor',[.7 .7 .7])
fill([Rt_epiestim_true(1).t_end'; flipud(Rt_epiestim_true(1).t_end'); Rt_epiestim_true(1).t_end(1)],[Rt_true_epiestim_ci(:,1,scenario); flipud(Rt_true_epiestim_ci(:,2,scenario)); Rt_true_epiestim_ci(1,1,scenario)],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
plot(Rt_epiestim_true(1).t_end',Rt_true_epiestim(:,scenario),'k','LineWidth',1)
plot(Rt_epiestim_obs(1).t_end',Rt_obs_epiestim(:,scenario),'--k','LineWidth',1)
plot((0:round(tmax*dt))',ones(round(tmax*dt)+1,1),'Color',[.6 .6 .6])
plot([20 20],[0 6],'k','LineWidth',.75)
ylabel('R_t')
if scenario==1
legend('True','Observed') 
end
xlim([20 70])
ylim([0 6])

subplot(5,4,4*scenario)
hold on
bar(Rt_epiestim_true(1).t_end(1)-1+find(bias_epiestim_lowerci(:,scenario)>0),bias_epiestim_lowerci(bias_epiestim_lowerci(:,scenario)>0,scenario),'k')
bar(Rt_epiestim_true(1).t_end(1)-1+find(bias_epiestim_upperci(:,scenario)<0),bias_epiestim_upperci(bias_epiestim_upperci(:,scenario)<0,scenario),'k')
xlim([20 70])
ylim([-2 2])
if scenario==1
legend('Bias')
end
end

%%
gtext('A)','FontWeight','bold','FontSize',8)
gtext('B)','FontWeight','bold','FontSize',8)
gtext('C)','FontWeight','bold','FontSize',8)
gtext('D)','FontWeight','bold','FontSize',8)
gtext('E)','FontWeight','bold','FontSize',8)
%gtext('F)','FontWeight','bold','FontSize',8)

%% Figure 2 main text

figure
for scenario=6:11
subplot(6,4,4*(scenario-6)+1)
%semilogy(cumsum(Npositive(:,scenario)),'b') 
%hold on
%semilogy(cumsum(Ntested(:,scenario)),'r')
%semilogy(cumsum(NewCase(:,scenario)),'--k')
hold on
plot(Ntested(:,scenario),'--k','LineWidth',1)
plot(Npositive(:,scenario),'k','LineWidth',1)
if scenario==6
    legend('Tests','Confirmed cases')
end
ylabel('# of cases')
xlim([20 70])
ylim([0 30000])

subplot(6,4,4*(scenario-6)+2)
plot(100*Ppositive(:,scenario),'k')
ylabel('% positive')
ylim([0 80])
xlim([20 70])
if scenario==6
    legend('%_{pos}')
end
if scenario==11
    xlabel('Days since start of epidemic')
end

subplot(6,4,4*(scenario-6)+3)
hold on
plot(Rt_epiestim_true(1).t_end',Rt_true_epiestim(:,scenario),'k','LineWidth',1)
plot(Rt_epiestim_obs(1).t_end',Rt_obs_epiestim(:,scenario),'--k','LineWidth',1)
fill([Rt_epiestim_obs(1).t_end'; flipud(Rt_epiestim_obs(1).t_end'); Rt_epiestim_obs(1).t_end(1)],[Rt_obs_epiestim_ci(:,1,scenario); flipud(Rt_obs_epiestim_ci(:,2,scenario)); Rt_obs_epiestim_ci(1,1,scenario)],[.7 .7 .7],'EdgeColor',[.7 .7 .7])
fill([Rt_epiestim_true(1).t_end'; flipud(Rt_epiestim_true(1).t_end'); Rt_epiestim_true(1).t_end(1)],[Rt_true_epiestim_ci(:,1,scenario); flipud(Rt_true_epiestim_ci(:,2,scenario)); Rt_true_epiestim_ci(1,1,scenario)],[.9 .9 .9],'EdgeColor',[.9 .9 .9])
plot(Rt_epiestim_true(1).t_end',Rt_true_epiestim(:,scenario),'k','LineWidth',1)
plot(Rt_epiestim_obs(1).t_end',Rt_obs_epiestim(:,scenario),'--k','LineWidth',1)
plot((0:round(tmax*dt))',ones(round(tmax*dt)+1,1),'Color',[.6 .6 .6])
plot([20 20],[0 8],'k','LineWidth',.75)
ylabel('R_t')
if scenario==6
legend('True','Observed') 
end
xlim([20 70])
ylim([0 8])

subplot(6,4,4*(scenario-5))
hold on
bar(Rt_epiestim_true(1).t_end(1)-1+find(bias_epiestim_lowerci(:,scenario)>0),bias_epiestim_lowerci(bias_epiestim_lowerci(:,scenario)>0,scenario),'k')
bar(Rt_epiestim_true(1).t_end(1)-1+find(bias_epiestim_upperci(:,scenario)<0),bias_epiestim_upperci(bias_epiestim_upperci(:,scenario)<0,scenario),'k')
xlim([20 70])
ylim([-5 5])
if scenario==6
legend('Bias')
end

end

%% Figure 3 main text (US data)
figure
subplot(1,4,1)
hold on
[ax,y1,y2]=plotyy(dateUS,newtests,dateUS,newcases);
set(y1,'Color',[.7 .7 .7],'LineWidth',1)
set(y2,'Color','k','LineWidth',1)
set(ax(1),'XLim',[dateUS(1)-3 dateUS(end)+2],'XTick',[dateUS(1)-3; dateUS([59 120],1); dateUS(end)+2],'XTickLabel',{'Mar-01','May-01','Jul-01','Sep-01'},'YColor','k','FontSize',8)
set(ax(2),'XLim',[dateUS(1)-3 dateUS(end)+2],'XTick',[dateUS(1)-3; dateUS([59 120],1); dateUS(end)+2],'XTickLabel',{'Mar-01','May-01','Jul-01','Sep-01'},'YColor','k','FontSize',8)
ax(1).YLabel.String='# of tests';
ax(2).YLabel.String='# of positive cases';
legend('Tests','Positive cases')
    
subplot(1,4,2)
%semilogy(dateUS(1:29),ntest(1:29),'r')
semilogy(dateUS(1:29),newtests(1:29),'Color',[.7 .7 .7],'LineWidth',1)
hold on
%semilogy(dateUS(1:29),npos(1:29),'b')
semilogy(dateUS(1:29),newcases(1:29),'k','LineWidth',1)
semilogy([dateUS(21) dateUS(21)],[1 1e7],':k')
semilogy(dateUS(1:29),glmval(usgrowthrate_test,1:29,'log'),'Color',[.7 .7 .7],'LineStyle','--')
semilogy(dateUS(1:29),glmval(usgrowthrate,1:29,'log'),'Color','k','LineStyle','--')
datetick('x','mmm-dd-yy')
ylabel('# of tests')
set(gca,'XTick',[dateUS(1)-3; dateUS([12 29],1)],'XTickLabel',{'Mar-01','Mar-15','Apr-01'},'FontSize',8)
xlim([dateUS(1)-3 dateUS(29)+1])
ylim([0 1e6])
legend('Tests','Positive cases','Location','NW')

subplot(1,4,3)
hold on
plot(dateUS,100*ppos,'Color',[.7 .7 .7])
plot(dateUS,100*ppos_mavg,'k','LineWidth',1)
ylabel('% positive')
legend('Daily data','15-day moving avg')
datetick('x','mmm-dd-yy')
xlim([dateUS(1)-3 dateUS(end)+2])
set(gca,'XTick',[dateUS(1)-3; dateUS([59 120],1); dateUS(end)+2],'XTickLabel',{'Mar-01','May-01','Jul-01','Sep-01'},'FontSize',8)

subplot(1,4,4)
hold on
plot(dateUS(1)-1+Rt_epiestim_US.t_end',Rt_US_epiestim,'k','LineWidth',1)
fill([dateUS(1)-1+Rt_epiestim_US.t_end'; dateUS(1)-1+flipud(Rt_epiestim_US.t_end'); dateUS(1)-1+Rt_epiestim_US.t_end(1)],[Rt_US_epiestim_ci(:,1); flipud(Rt_US_epiestim_ci(:,2)); Rt_US_epiestim_ci(1,1)],[.7 .7 .7],'EdgeColor',[.7 .7 .7])
plot(dateUS(1)-1+Rt_epiestim_US.t_end',Rt_US_epiestim,'k','LineWidth',1)
%plot(dateUS(1:end-12,1),Rt_US(1:end-11,1),'Color',[0 .7 0],'LineWidth',1)
%plot(dateUS(1:end-1,1),Rt_US(1:end,1),':k')
%plot(dateUS(1:end-1,1),Rt_US_adj(1:end,1),'k')
%fill([dateUS(1:end-1,1); dateUS(end-1:-1:1,1); dateUS(1)],[Rt_US_ci(1:end,1); Rt_US_ci(end:-1:1,2); Rt_US_ci(1,1)],[.7 1 .7],'EdgeColor',[.7 1 .7])
%plot(dateUS(1:end-1,1),Rt_US(1:end,1),':k')
%plot(dateUS(1:end-1,1),Rt_US_adj(1:end,1),'k')
plot(dateUS,ones(length(dateUS),1),'--k')
%plot(dateUS(1:end-12,1),Rt_US(1:end-11,1),'Color',[0 .7 0],'LineWidth',1)
ylabel('R_t')
datetick('x','mmm-dd-yy')
set(gca,'XTick',[dateUS(1)-3; dateUS([59 120],1); dateUS(end)+2],'XTickLabel',{'Mar-01','May-01','Jul-01','Sep-01'},'FontSize',8)
xlim([dateUS(1)-3 dateUS(end)+2])
ylim([0 5])
legend('Estimated \it{R_t}','95% CI') %,'Uncorrected','Corrected')
%%
gtext('A)','FontWeight','bold','FontSize',8)
gtext('B)','FontWeight','bold','FontSize',8)
gtext('C)','FontWeight','bold','FontSize',8)
gtext('D)','FontWeight','bold','FontSize',8)

%% Figure 4 (final figure for select states)
%Divide into early vs late periods
st=[5 4 19 10 23 35];
n=length(st);
early=29;

figure
for i=1:n
    subplot(4,n,i)
    %semilogy(dateUS,npos_bystate(:,st(i)),'b')
    %hold on
    %semilogy(dateUS,ntest_bystate(:,st(i)),'r')
    [ax,y1,y2]=plotyy(dateUS(1:early,1),newtests_bystate(1:early,st(i))/1000,dateUS(1:early,1),newcases_bystate(1:early,st(i))/1000);
    set(y1,'Color',[.7 .7 .7],'LineWidth',1)
    set(y2,'Color','k','LineWidth',1)
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
    [ax,y1,y2]=plotyy(dateUS(1)-1+Rt_epiestim_bystate(st(i)).t_end',Rt_bystate_epiestim(:,st(i)),dateUS(1:early,1),100*ppos_bystate(1:early,st(i)));
    set(y1(1),'Color','k','LineWidth',1)
    set(y2(1),'Color','k','LineStyle','--')
    set(ax(1),'XLim',[dateUS(1)-3 dateUS(early)+1],'XTick',[dateUS(1)-3; dateUS(29,1)],'XTickLabel',{'Mar-01','Apr-01'},'YLim',[0 10],'YTick',0:5:10,'YColor','k','FontSize',8)
    set(ax(2),'XLim',[dateUS(1)-3 dateUS(early)+1],'XTick',[dateUS(1)-3; dateUS(29,1)],'XTickLabel',{'Mar-01','Apr-01'},'YLim',[0 100],'YTick',0:50:100,'YColor','k','FontSize',8)
    if i==1
        ax(1).YLabel.String='R_t';
    elseif i==n    
        ax(2).YLabel.String='% positive';
        legend('Estimated R_t','Percent positive')
    end
    fill([dateUS(1)-1+Rt_epiestim_bystate(st(i)).t_end'; dateUS(1)-1+flipud(Rt_epiestim_bystate(st(i)).t_end'); dateUS(1)-1+Rt_epiestim_bystate(st(i)).t_end(1)],[Rt_bystate_epiestim_ci(:,1,st(i)); flipud(Rt_bystate_epiestim_ci(:,2,st(i))); Rt_bystate_epiestim_ci(1,1,st(i))],[.7 .7 .7],'EdgeColor',[.7 .7 .7])
    plot(dateUS(1)-1+Rt_epiestim_bystate(st(i)).t_end',Rt_bystate_epiestim(:,st(i)),'k','Linewidth',1)
    plot(dateUS(1:early,1),ones(early,1),'Color',[.7 .7 .7])
    
    subplot(4,n,2*n+i)
    %semilogy(dateUS,npos_bystate(:,st(i)),'b')
    %hold on
    %semilogy(dateUS,ntest_bystate(:,st(i)),'r')
    [ax,y1,y2]=plotyy(dateUS(early+1:end,1),newtests_bystate(early+1:end,st(i))/1000,dateUS(early+1:end,1),newcases_bystate(early+1:end,st(i))/1000);
    set(y1,'Color',[.7 .7 .7],'LineWidth',1)
    set(y2,'Color','k','LineWidth',1)
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
    [ax,y1,y2]=plotyy(dateUS(1)-1+Rt_epiestim_bystate(st(i)).t_end',Rt_bystate_epiestim(:,st(i)),dateUS(early+1:end,1),100*ppos_bystate(early+1:end,st(i)));
    set(y1(1),'Color','k','LineWidth',1)
    set(y2(1),'Color','k','LineStyle','--')
    set(ax(1),'XLim',[dateUS(early) dateUS(end)+1],'XTick',dateUS([59 120],1),'XTickLabel',{'May-01','Jul-01'},'YLim',[0 2],'YTick',0:2,'YColor','k','FontSize',8)
    set(ax(2),'XLim',[dateUS(early) dateUS(end)+1],'XTick',dateUS([59 120],1),'XTickLabel',{'May-01','Jul-01'},'YLim',[0 40],'YTick',0:20:40,'YColor','k','FontSize',8)
    if i==1
        ax(1).YLabel.String='R_t';
    elseif i==n    
        ax(2).YLabel.String='% positive';
    end
    fill([dateUS(1)-1+Rt_epiestim_bystate(st(i)).t_end'; dateUS(1)-1+flipud(Rt_epiestim_bystate(st(i)).t_end'); dateUS(1)-1+Rt_epiestim_bystate(st(i)).t_end(1)],[Rt_bystate_epiestim_ci(:,1,st(i)); flipud(Rt_bystate_epiestim_ci(:,2,st(i))); Rt_bystate_epiestim_ci(1,1,st(i))],[.7 .7 .7],'EdgeColor',[.7 .7 .7])
    plot(dateUS(1)-1+Rt_epiestim_bystate(st(i)).t_end',Rt_bystate_epiestim(:,st(i)),'k','Linewidth',1)
    plot(dateUS(early+1:end),ones(length(dateUS)-early,1),'Color',[.7 .7 .7])
    plot([dateUS(early) dateUS(early)],[0 2],'k')
    
end

%%
gtext('A)','FontWeight','bold','FontSize',8)
gtext('B)','FontWeight','bold','FontSize',8)
gtext('C)','FontWeight','bold','FontSize',8)
gtext('D)','FontWeight','bold','FontSize',8)
