%function Calc_MMBF
% Author: Shelley Zhang
% Email:  szhang@ltk.com
% Date:   Oct 2017
% Purpose:Phase II report, based on daily mileage data, failures, and work order counts
%% To load Miles and Incident data
clear all
close all
folder = 'P:\C4314 MBTA-X29PS01 Evaluation of Light Rail Vehicles\Type 8 Reliability Analysis\MCRS\DATA_mat';
folder = 'C:\Users\SZhang\Documents\MATLAB\T8\DATA_mat';
date_first = '2015-07-01';
date_last = '2018-12-31'; 
date_range = [date_first '_to_' date_last];
file_inc = ls([folder '\Incidents_' date_range '.mat'])
load([folder '\' file_inc])
file_mileage = ls([folder '\Mileages_' date_range '.mat'])
load([folder '\' file_mileage])
units = length(Mileage.Unit);
fprintf('Completed Reading in Incident and Mileage data\n');
Invertals = size(Mileage.Vals,2);
%Mileage.Vals = [Mileage.Vals zeros(units,ceil(Invertals/3)*3-Invertals)];


%% To bin into quarters and select interrupted incidents only
% may need further exclude Z98 accidents or Z20 Dirty Interior in Incident.Content(:,5)
%dates = char(zeros(length(Incident.Content),20));
dates_inc = datestr(Incident.Content(:,3));
dates_inc = datetime(dates_inc);
incident_in_quarter = my_discretize(datetime(datestr(Mileage.start_date)),dates_inc,'quarter');    
category_in_quarter = discretize(dates_inc,'quarter','categorical'); %output Q2 2010, Q3 2010 etc

incident_in_month = my_discretize(datetime(datestr(Mileage.start_date)),dates_inc,'month');    

months_inc = max(incident_in_month);
quarters_inc = max(incident_in_quarter);
years_inc = ceil(quarters_inc/4);

quarterfiller  = zeros(years_inc*4-quarters_inc, 1)'; % in case quarters are not integers of 4
monthfiller = zeros(quarters_inc*3-months_inc,1)';
%label_year = [repmat('Year ',[years 1]) num2str(2009+[1:years]')];
%Fiscal year 2011 is Jul 2011 -Jun 2012
label_year = mat2cell((year(datetime(datestr(Mileage.start_date)))+[1:years_inc]),1,ones(1,years_inc));
label_year{end}=[num2str(label_year{end}) ' to date'];
label_quarter = char(unique(category_in_quarter));
%Trip Lost Service Interruption at 14th column
incident_intrpt_index = Incident.events_intrpt;
%to include all incidents
%incident_intrpt_index = true(size(incident_in_quarter));
incident_intrpt_in_quarter = incident_in_quarter(incident_intrpt_index);
h = histogram(incident_intrpt_in_quarter,0.5:(quarters_inc+0.5)); 
fleet_incident_intrpt_counts_quarter = h.Values;
close
%if counting all incidents
h2 = histogram(incident_in_quarter,0.5:(quarters_inc+0.5));
fleet_incident_counts_quarter = h2.Values;
close
%% Percentage of failures that result in Removal from Service by Year
% incident_rfs_index = (strcmp(Incident.Content(:,7),'UNLOAD/REMOVE FROM SERVCE'));
% incident_rfs_in_quarter = incident_in_quarter(incident_rfs_index);
% h = histcounts(incident_rfs_in_quarter,0.5:4:(quarters_inc+0.5));
% h2 = histcounts(incident_in_quarter,0.5:4:(quarters_inc+0.5));
% figure,bar([h2-h;h]','stacked'),colormap(specialcmap('gkb'))
% text(1:years_inc, h2, num2str((100*h./h2)','%2.1f percent'),'HorizontalAlignment','center','VerticalAlignment','bottom')
% xticks(1:years_inc)
% set(gca,'XTickLabel',label_year)
% xlabel('Intervals [Year]')
% ylabel('Incident [Count]')
% title('Incidents that Result in UNLOAD/REMOVE FROM SERVICE by Year')

%% Incident Counts by car
%Car Unit in Incident.Content(:,1)
incident_intrpt_by_car = str2double(Incident.Content(incident_intrpt_index,1));
cars_with_intrpt_incidents = unique(incident_intrpt_by_car);
cars_with_incidents = unique(str2double(Incident.Content(:,1)));
validcars = ismember(Mileage.Unit,cars_with_intrpt_incidents);
validcars2 = ismember(Mileage.Unit,cars_with_incidents);

%%%%%%%%%%%%%Out of Service Cars%%%%%%%%%%%%%%%%
cars_oos = [3803 3807 3808 3832 3836 3854 3873 3879 3891 3893]; %3879 not even show up in Inspections Due
% validcars(ismember(Mileage.Unit,cars_oos))=0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cars_without_intrpt_incidents = Mileage.Unit(setdiff([1:units],find(validcars)));
%initialize
incident_intrpt_by_car_by_quarter = zeros(units,quarters_inc);
% form a bivariate histogram of carunit by quarter
incident_intrpt_by_car_by_quarter(validcars,:) = histcounts2(incident_intrpt_by_car, incident_intrpt_in_quarter,...
    [cars_with_intrpt_incidents-0.5;(cars_with_intrpt_incidents(end)+0.5)], [0.5:1:(quarters_inc+0.5)] );

%initialize
incident_intrpt_by_car_by_year = zeros(units,years_inc);
% form a bivariate histogram of carunit by year
incident_intrpt_by_car_by_year(validcars,:) = histcounts2(incident_intrpt_by_car, incident_intrpt_in_quarter,...
    [cars_with_intrpt_incidents-0.5;(cars_with_intrpt_incidents(end)+0.5)], [0.5:4:(years_inc*4+0.5)]  );

% to count overall incidents by car
incident_count_by_car = zeros(units,1);
incident_count_by_car(validcars,:) = histcounts(incident_intrpt_by_car,...
    [cars_with_intrpt_incidents'-0.5 max(cars_with_intrpt_incidents)+0.5]);
%initialize
incident_count_by_car_by_quarter = zeros(units,quarters_inc);
incident_count_by_car_by_quarter(validcars2,:) = histcounts2(str2double(Incident.Content(:,1)), incident_in_quarter,...
    [cars_with_incidents-0.5;(cars_with_incidents(end)+0.5)], [0.5:1:(quarters_inc+0.5)] );


%% To calculate MMBF for the whole fleet by Quarters
% quarterly arrangement
firstday_in_quarter = datenum(datetime(datestr(Incident.start_date))+calquarters(0:quarters_inc));

%firstday_in_quarter = datenum(0+calmonths(0:3:years_inc*12));
miles_by_car_by_quarter = zeros(units,quarters_inc);
%for q=1:quarters_inc
    %miles_by_car_by_quarter(:,q) = sum(Mileage.Vals(:,(firstday_in_quarter(q)+1:firstday_in_quarter(q+1))-Incident.start_date),2);
%end
miles_by_car_by_quarter = squeeze(sum(reshape([Mileage.Vals repmat(monthfiller,units,1)],units,3,quarters_inc),2));
fleet_miles = sum(miles_by_car_by_quarter,1);
MMBF_quarterly = fleet_miles./fleet_incident_intrpt_counts_quarter;
MMBI_quarterly = fleet_miles./fleet_incident_counts_quarter;
figure,bar(MMBF_quarterly,'FaceColor',[0 .5 .5],'EdgeColor',[0 .5 .5],'LineWidth',1.5)
yneg = MMBF_quarterly-fleet_miles./(fleet_incident_intrpt_counts_quarter+2);
ypos = fleet_miles./(fleet_incident_intrpt_counts_quarter-2)-MMBF_quarterly;
hold on, errorbar(1:quarters_inc, MMBF_quarterly,yneg,ypos,'r.','Linewidth',1.5);
xticks(1:quarters_inc)
set(gca,'XTickLabel',label_quarter)
set(gca,'XTickLabelRotation',90)
xlabel('Intervals [Quarter]')
ylabel('MMBF [Mile]')
title('Fleet Mean Miles Between Failures by Quarter')
text(1:quarters_inc,MMBF_quarterly+250,num2str(MMBF_quarterly', '%.f'),'HorizontalAlignment','center')
MMBF_quarterly'
pause
MMBI_quarterly'
pause
[Bm,Im]=sort(miles_by_car_by_quarter(:,end),'descend');
[Mileage.Unit(Im)' Bm]

%% yearly arrangement
%%%%Be careful of where quarterfill is placed in front of at back of fleet_miles
MMBF_yearly = sum(reshape([fleet_miles quarterfiller],4,years_inc),1)./sum(reshape([fleet_incident_intrpt_counts_quarter quarterfiller],4,years_inc),1);
MMBI_yearly = sum(reshape([fleet_miles quarterfiller],4,years_inc),1)./sum(reshape([fleet_incident_counts_quarter quarterfiller ],4,years_inc),1);
% MMBF_yearly = sum(reshape([fleet_miles(3:end)],4,years_inc),1)./sum(reshape([fleet_incident_counts_quarter(3:end)],4,years_inc),1);
figure, bar(MMBF_yearly,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)
xticks(1:years_inc)
set(gca,'XTickLabel',label_year)
xlabel('Intervals [Year]')
ylabel('MMBF [Mile]')
title('Fleet Mean Miles Between Failures by Year')
text(1:years_inc, MMBF_yearly'-500,num2str(MMBF_yearly','%0.0f'),'color','white',...
    'HorizontalAlignment','center','VerticalAlignment','bottom')
%For executive summary
MMBF_allyears = sum([quarterfiller fleet_miles],2)/sum([quarterfiller fleet_incident_counts_quarter],2);
MMBF_last3years = sum(fleet_miles(end-11:end),2)/sum(fleet_incident_counts_quarter(end-11:end),2);
% fprintf('MMBF all years: %.1f, MMBF last 3 years: %.1f\n',MMBF_allyears,MMBF_last3years)
% hold on, plot(0.5:years_inc+0.5,MMBF_allyears*ones(1,years_inc+1),'k','linewidth',6)
% plot(years_inc-2.5:years_inc+0.5,MMBF_last3years*ones(1,4),'b','linewidth',6)
% plot(years_inc-0.5:0.1:years_inc+0.5,MMBF_yearly(end)*ones(1,11),'g','linewidth',6)
% text((years_inc+1)*ones(1,3), [MMBF_allyears MMBF_last3years MMBF_yearly(end)],[num2str([MMBF_allyears MMBF_last3years MMBF_yearly(end)]','%0.0f') repmat(' Miles',[3 1])],'HorizontalAlignment','center','VerticalAlignment','bottom')
% set(gca,'box','off')
%% MMBF per car over all time, excludes non-service car, i.e. travel miles less than 50 quarterly
%car_incident_counts(car_incident_counts<40)=0;
MMBF_car = sum(Mileage.Vals,2)./incident_count_by_car;
%MMBF_car(ismember(Mileage.Unit,cars_oos))=Inf; 
[MMBF_car_sorted,I]=sort(MMBF_car);
infc = length(find(isinf(MMBF_car_sorted)));
figure, barh(MMBF_car_sorted,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)
yticks(1:units-infc),ylim([0 units-infc+1])
set(gca,'YTickLabel',Mileage.Unit(I(1:end-infc)))
%set(gca,'YTickLabelRotation',90)
ylabel('Car [Unit]')
xlabel('MMBF [Mile]')
title(['Individual Car Mean Miles Between Failures over ' datestr(Mileage.start_date) ' to ' datestr(Mileage.end_date)])
plotmeanstd(MMBF_car_sorted,units-infc)
carrank_MMBF_alltime = setdiff(Mileage.Unit(I(1:end-infc))',cars_oos','stable')
%% MMBF by car by quarter
close all
%%%%%%%%%%%% assuming less than 50 miles within a quarter means that the
%%%%%%%%%%%% car is not in service
%miles_by_car_by_quarter(miles_by_car_by_quarter<50)=0;
%%%%%%%%%%%%%%%%%%%%
% MMBF_car_quarter = miles_car_quarter(validcars,:)./incident_interrupted_car_quarter;
%incident_interrupted_car_quarter(incident_interrupted_car_quarter==0)=1;
%incident_intrpt_by_car_by_quarter(incident_intrpt_by_car_by_quarter==0)=1;
tmp = incident_intrpt_by_car_by_quarter;
tmp(tmp==0)=1; 
MMBF_car_quarter = miles_by_car_by_quarter./tmp;
[B,I]=sort(MMBF_car_quarter(:,end),'descend');
[Mileage.Unit(I)' miles_by_car_by_quarter(I,end) incident_intrpt_by_car_by_quarter(I,end) B]
pause
clear tmp
tmp = incident_count_by_car_by_quarter;
tmp(tmp==0)=1; 
MMBI_car_quarter = miles_by_car_by_quarter./tmp;
[B,I]=sort(MMBI_car_quarter(:,end),'descend');
[Mileage.Unit(I)' miles_by_car_by_quarter(I,end) incident_count_by_car_by_quarter(I,end) B]

pause
clear tmp
tmp = sum(incident_intrpt_by_car_by_quarter(:,end-4:end-1),2);
tmp(tmp==0)=1;
MMBF_car_2018 = sum(miles_by_car_by_quarter(:,end-4:end),2)./tmp;
[B,I]=sort(MMBF_car_2018);
[Mileage.Unit(I)' sum(miles_by_car_by_quarter(I,end-4:end-1),2) tmp(I) B]


pause
clear tmp
tmp = sum(incident_count_by_car_by_quarter(:,end-4:end-1),2);
tmp(tmp==0)=1;
MMBI_car_2018 = sum(miles_by_car_by_quarter(:,end-4:end),2)./tmp;
[B,I]=sort(MMBI_car_2018);
[Mileage.Unit(I)' sum(miles_by_car_by_quarter(I,end-4:end-1),2) tmp(I) B]
pause
%MMBF_car_quarter(isinf(MMBF_car_quarter)|isnan(MMBF_car_quarter))=0;
maxMMBFq = max(MMBF_car_quarter(MMBF_car_quarter<Inf));
figure,imagesc(miles_by_car_by_quarter),colormap(specialcmap('summer')),title('Mileage by Quarter (-X) by Car(-Y)')
ylabel('Car [Unit]')
yticks([1:5:units units])
tmp = get(gca,'YTick')+Mileage.Unit(1)-1; 
set(gca,'YTickLabel',tmp,'fontsize',18,'fontweight','bold')
xticks(1:quarters_inc),set(gca,'XTickLabel',label_quarter,'XTickLabelRotation',90)
colorbar

% figure,imagesc(MMBF_car_quarter),colormap(specialcmap('okg',6850/max(MMBF_car_quarter(:)))),title('MMBF by Quarter (-X) by Car(-Y)')
% ax=gca;
% ylabel('Car [Unit]')
% yticks([1:units])%yticks([1:5:units units])
% tmp = get(gca,'YTick')+Mileage.Unit(1)-1; 
% set(gca,'YTickLabel',tmp,'fontsize',18,'fontweight','bold')
% yrule=ax.YAxis;yrule.FontSize=9;yrule.FontWeight='bold';yrule.TickLabelRotation=0;
% xticks(1:quarters_inc),set(gca,'XTickLabel',label_quarter,'XTickLabelRotation',15)
% colorbar
% 
% 
% figure,imagesc(incident_intrpt_by_car_by_quarter),colormap(specialcmap('summer')),title('Incident Counts by Quarter (-X) by Car (-Y)')
% ylabel('Car [Unit]')
% yticks([1:5:units units])
% tmp = get(gca,'YTick')+Mileage.Unit(1)-1; 
% set(gca,'YTickLabel',tmp,'fontsize',18,'fontweight','bold')
% xticks(1:quarters_inc),set(gca,'XTickLabel',label_quarter,'XTickLabelRotation',90)
% 
% colorbar

%% Plot MMBF by Quarter for Every Unit
screensize=get(0,'ScreenSize');
figure('Position',screensize)
for u = 1:units %Weekly Mileage
    %plot(MMBF_car_quarter(u,:),'k.-','MarkerSize',40,'Linewidth',3,'Color',[0 .5 .5])
    bar(MMBF_car_quarter(u,:),'FaceColor',[0 .5 .5],'FaceAlpha',0.5)
    zero_inc = find(incident_intrpt_by_car_by_quarter(u,:)==0);
    if ~isempty(zero_inc)
        hold on, plot(zero_inc,miles_by_car_by_quarter(u,zero_inc),'O','MarkerSize',20,'Linewidth',3,'Color',[.2 .8 .2])
        text(zero_inc,MMBF_car_quarter(u,zero_inc)+350,'Inc=0','FontSize',16,'Rotation',90,'Color',[.2 .8 .2])
    end
    zero_miles = find(miles_by_car_by_quarter(u,:)==0);
    if ~isempty(zero_miles)
        text(zero_miles,MMBF_car_quarter(u,zero_miles)+300,'& Mile=0','FontSize',16,'Rotation',90)
    end
    hold on, 
    tmp = MMBF_car_quarter; tmp(tmp<50)=NaN;
    fleetmean = nanmean(tmp,1);
    plot(fleetmean,'.-','MarkerSize',40,'Linewidth',3,'Color',[0 0 0])
    %bar(fleetmean,'FaceColor',[0 .5 .5],'FaceAlpha',0.2)
    text(quarters_inc+0.1,fleetmean(end),'Fleet Average','FontSize',12,'FontWeight','bold')
    plot(fleetmean+1*std(MMBF_car_quarter,1),'.--','MarkerSize',40,'Linewidth',3,'Color',[0 0 0])
    text(quarters_inc+0.1,fleetmean(end)+std(MMBF_car_quarter(:,end),1),'Plus 1SD','FontSize',12,'FontWeight','bold')
    plot(fleetmean-1*std(MMBF_car_quarter,1),'.--','MarkerSize',40,'Linewidth',3,'Color',[0 0 0])
    text(quarters_inc+0.1,fleetmean(end)-std(MMBF_car_quarter(:,end),1),'Minus 1SD','FontSize',12,'FontWeight','bold')
    xticks(1:quarters_inc);
    set(gca,'XTickLabel',label_quarter,'XTickLabelRotation',17,'FontSize',20,'FontWeight','bold')
    %ax=gca;xrule=ax.XAxis;xrule.FontSize=16;xrule.FontWeight='bold';xrule.TickLabelRotation=0;
    ylabel('MMBF [Mile]','FontSize',20,'FontWeight','bold')
    title(['Unit 0' num2str(Mileage.Unit(u)) ' Quarterly MMBF Bargraph over Last Three Years'])
    xlim([0 quarters_inc+1])
    ylim([0 ceil(maxMMBFq/500)*500])
    mkdir([folder '..\..\figs\UnitFigs\',num2str(Mileage.Unit(u))])
    saveas(gcf,[[folder '..\..\figs\UnitFigs\',num2str(Mileage.Unit(u))] '\' num2str(Mileage.Unit(u)) '_QuarterlyMMBF'...
        datestr(Mileage.start_date) 'to' datestr(Mileage.end_date) '.png'])
    hold off
    
end
%% MMBF by car by year
miles_car_year = squeeze(sum(reshape([miles_by_car_by_quarter repmat(quarterfiller,[units 1])],units,4,years_inc),2));
miles_car_year(miles_car_year<200) = 0;
MMBF_car_year = miles_car_year./incident_intrpt_by_car_by_year;
% To clear NaN and Inf as 0/0 = NaN, number/0 = Inf;
MMBF_car_year(isinf(MMBF_car_year)|isnan(MMBF_car_year))=0;
%for exec summary 
% MMBF_car_last3years = sum(miles_car_year(validcars,end-2:end),2)./sum(incident_intrpt_by_car_by_year(:,end-2:end),2);
% [sorted3,sind3y] = sort(MMBF_car_last3years);

% MMBF_car_last3years = sum(miles_car_year(validcars,end-2:end),2)./sum(incident_intrpt_by_car_by_year(validcars,end-2:end),2);
% [sorted3,sind3y] = sort(MMBF_car_last3years);
% 
% MMBF_car_last7years = sum(miles_car_year(validcars,:),2)./sum(incident_intrpt_by_car_by_year(validcars,end-6:end),2);
% [sorted7,sind7y] = sort(MMBF_car_last7years);

%lastyears = input(['Please input the number of latest years to plot (1-' num2str(years)],'s');
lastyears = 2
toplot = miles_car_year(validcars,:);
toplot = incident_intrpt_by_car_by_year;
toplot = MMBF_car_year(validcars,:);
figure, 
for n=1:lastyears
    tmp = (lastyears-n);
    subplot(1,lastyears,n),barh(Mileage.Unit(validcars),toplot(:,end-tmp),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5),
    title(label_year(end-tmp),'Fontsize',16)
    %ylim([min(nonzeros(toplot)) max(nonzeros(toplot))])
    set(gca,'XLim',[0 max(toplot(:))])
    xlabel('MMBF [Mile]','FontSize',14,'fontweight','bold')
    yticks(cars_with_intrpt_incidents)
    yticklabels(Mileage.Unit(validcars))
    ylim([min(cars_with_intrpt_incidents)-1 max(cars_with_intrpt_incidents)+1])
    %axis tight
end
figure %to sort it
for n=1:lastyears
    tmp = (lastyears-n);
    clear sorted*
    [sorted,sortedi]=sort(toplot(:,end-tmp),1,'descend');
    tmp2 = Mileage.Unit(validcars);
    subplot(1,lastyears,n),barh(sorted,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5),
    title(label_year(end-tmp),'Fontsize',16)
    %ylim([min(nonzeros(toplot)) max(nonzeros(toplot))])
    set(gca,'XLim',[0 max(toplot(:))])
    xlabel('MMBF [Mile]','FontSize',14,'fontweight','bold')
    yticks(1:length(cars_with_intrpt_incidents))
    yticklabels(tmp2(sortedi))
    ylim([0 length(cars_with_intrpt_incidents)+1])
    %axis tight
end
carrank_MMBF_last_year=setdiff(tmp2(sortedi)',cars_oos','stable');

