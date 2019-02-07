function ExportVH_in_system
clear all
close all

folder = 'C:\Users\SZhang\Documents\MATLAB\T8\DATA_mat\';
filelist = ls([folder 'VehicleHistory*.mat']);
load([folder  filelist])
%clear out all work orders under Y- category (i.e. inspection related)
Maintenance.Content(startsWith(Maintenance.Content(:,9),'Y-'),:)=[];
%% Select only between desired start and end dates
desired_start_date = max(datetime(2015,7,1),datetime(datestr(Maintenance.start_date)));
desired_end_date = min(datetime(2018,8,15),datetime(datestr(Maintenance.end_date)));
% months_tot = calmonths(between(desired_start_date-1,desired_end_date));
dates_inc = datetime(Maintenance.Content(:,1));

workorder_in_month = my_discretize(datetime(datestr(Maintenance.start_date)),dates_inc,'month');    
months_inc = max(workorder_in_month);

tf=isbetween(Maintenance.Content(:,1),desired_start_date,desired_end_date);
%uint32 cuts memory size to 50% compared to double. uint32 stores 0 to 2^32-1=4e9 
tfi = uint32(find(tf==1)); 
event_symptom_code=char(Maintenance.Content(tf,9));

systems = unique(event_symptom_code(:,1));
opt = input(['System to request: D/[' systems' ']:'],'s');
if isempty(opt)
    opt = 'D';
end
if length(tfi)< 2^16
    rows1 = uint16(strfind(event_symptom_code(:,1)',opt));
    %rows where only letter space ( show up, i.e. high level system codes
    rows2 = uint16(strfind(reshape(event_symptom_code(:,2:3)',[1 length(event_symptom_code)*2]),'-:')+1)/2;
else
    rows1 = uint32(strfind(event_symptom_code(:,1)',opt));
    %rows where only letter space ( show up, i.e. high level system codes
    rows2 = uint32(strfind(reshape(event_symptom_code(:,2:3)',[1 length(event_symptom_code)*2]),'-:')+1)/2;
end
system_label = unique(event_symptom_code(rows2,:),'rows');
tmp = char(system_label);
system_opt = strtrim(replace(system_label(strfind(tmp(:,1)',opt),:),':',''));
filename = [system_opt '_Maintenance_' datestr(desired_start_date,'yyyy-mm-dd') '_to_' datestr(desired_end_date,'yyyy-mm-dd') '.xlsx']
inputcon = Maintenance.Content(tfi(rows1),:);
xlswrite([folder filename],Maintenance.Fieldname)
xlswrite([folder filename],inputcon,1,'A2')
%%
%clear Main* 
[tasks_unique,ia,ic]= unique(join([inputcon(:,4) inputcon(:,9)]),'stable');

tasks_laborhours = accumarray(ic,cell2mat(inputcon(:,5)));
%no tasks haing 0 hours
% %remove tasks that have 0 hours
% tasks_null = find(tasks_laborhours==0);
% tasks_laborhours(tasks_null)=[];
% %remove index of tasks that have 0 hours
% ia(tasks_null)=[];

tasks_dates_num = datenum(inputcon(ia,1));
tasks_symptoms = categorical(inputcon(ia,9));
tasks_symptoms_num = double(tasks_symptoms);
symptom_labels = categories(categorical(inputcon(ia,9)));
xbin = datenum(desired_start_date+calmonths(0:1:months_inc));
ybin = double(1:max(tasks_symptoms_num)+1);
%% 2D histogram by time and symotoms
%histogram bin edges is left inclusive
%2D histogram to split task counts by month and symptom codes


h2 = histogram2(tasks_dates_num,tasks_symptoms_num,xbin,ybin)
title('Truck Tasks Counts by Month and by Symptoms')
set(gca,'XTIck',xbin(1:end-1),'XTickLabel', datestr(xbin(1:end-1)),'XTickLabelRotation',-20)
set(gca,'YTIck',ybin(1:end-1)+0.5,'YTickLabel',symptom_labels ,'YTickLabelRotation',60)

% first_half = 1:round((length(ybin)-1)/2);
% figure,h2=plot(h.Values(:,first_half),'O-'),title('Truck Task Counts by Month')
% hold on, plot(h.Values(:,max(first_half)+1:end),'O-.'),legend(symptom_labels)
% figure,h2=plot(h.Values,'.-','MarkerSize',14),title('Truck Task Counts by Month')
% ii = (1:20)';
% %tmp =jet(128);set(h2, {'color'}, num2cell(tmp(1:5:99,:), 2));
% set(h2, {'color'}, num2cell(jet(length(ii)), 2));
% set(gca,'XTIck',1:months_tot,'XTickLabel', datestr(xbin(1:end-1),'yyyy-mmm'),'XTickLabelRotation',45)
% [ta,tb]=max(h.Values,[],1);
% text(tb,ta,char(symptom_labels),'Rotation',0)
figure
plot(sum(h2.Values,2),'.-','MarkerSize',20,'LineWidth',3);%(bar(sum(h.Values,2)),
set(gca,'XTIck',1:months_inc,'XTickLabel', datestr(xbin(1:end-1),'yy-mmm'),'XTickLabelRotation',90)
%set(gca,'XTIck',1:3:months_tot,'XTickLabel', datestr(xbin(1:3:end-1),'yy-mmm'),'XTickLabelRotation',90)
title('Truck Task Counts by Month')
figure,
for n=1:length(ybin)-1
    subplot(4,ceil(length(ybin)/4),n),plot(h2.Values(:,n),'.-','MarkerSize',20,'LineWidth',3),title(symptom_labels{n})
    set(gca,'XTIck',1:3:months_inc,'XTickLabel', datestr(xbin(1:3:end-1),'yy-mmm'),'XTickLabelRotation',90)
    ymax = get(gca,'ylim');
    rectangle('Position',[1 0 12 ymax(2)],'FaceColor',[0 0 0 0.3])
end
%figure,imagesc(h.Values)

%% findgroups and splitapply by time and symotoms

tasks_dates_month = double(discretize(datetime(inputcon(ia,1)),'month','categorical'));
[G,ID1,ID2] = findgroups(tasks_dates_month,tasks_symptoms_num);
tasks_counts_by_month_by_symptom = zeros(max(ID1),max(ID2));
tasks_laborhours_by_month_by_symptom = zeros(max(ID1),max(ID2));
newsub = sub2ind([max(ID1) max(ID2)],ID1,ID2);

tmp = splitapply(@sum,tasks_laborhours,G);
tasks_laborhours_by_month_by_symptom(newsub) = tmp;

tmp = splitapply(@sum,ones(size(G)),G);
tasks_counts_by_month_by_symptom(newsub) = tmp;
hours_per_task=16;
yaxismax = 1.05*max([sum(tasks_counts_by_month_by_symptom,2)*hours_per_task;sum(tasks_laborhours_by_month_by_symptom,2)]);

figure
yyaxis left
plot(sum(tasks_counts_by_month_by_symptom,2),'.-','MarkerSize',28,'LineWidth',3);%(bar(sum(h.Values,2)),
set(gca,'XTIck',1:months_tot,'XTickLabel', datestr(xbin(1:end-1),'yy-mmm'),'XTickLabelRotation',90)
%set(gca,'XTIck',1:3:months_tot,'XTickLabel', datestr(xbin(1:3:end-1),'yy-mmm'),'XTickLabelRotation',90)
ylim([0 yaxismax/hours_per_task])
yyaxis right
plot(sum(tasks_laborhours_by_month_by_symptom,2),'.-','MarkerSize',28,'LineWidth',3,'Color',[0.85 0.33 0.1 0.5]);%(bar(sum(h.Values,2)),
set(gca,'XTIck',1:months_tot,'XTickLabel', datestr(xbin(1:end-1),'yy-mmm'),'XTickLabelRotation',90)
%set(gca,'XTIck',1:3:months_tot,'XTickLabel', datestr(xbin(1:3:end-1),'yy-mmm'),'XTickLabelRotation',90)
ylim([0 yaxismax])
ymax = get(gca,'ylim');
rectangle('Position',[1 0 14 ymax(2)],'FaceColor',[0 0 0 0.3])
title('Truck Task Counts (B) and Labor Hours (O) by Month')

figure,annotation('textbox',[0.4 1 0.2 0.05],'String','Task Counts (B) and Labor Hours (O) by Month',...
    'FitBoxToText','on','LineStyle','none');

for n=1:length(ybin)-1
    yaxismax = 1.1*max([tasks_counts_by_month_by_symptom(:,n)*hours_per_task;tasks_laborhours_by_month_by_symptom(:,n)]);
    subplot(5,4,n),title([symptom_labels{n}])
    yyaxis left
    plot(tasks_counts_by_month_by_symptom(:,n),'.-','MarkerSize',20,'LineWidth',3)
    ylim([0 yaxismax/hours_per_task])
    yyaxis right
    plot(tasks_laborhours_by_month_by_symptom(:,n),'.-','MarkerSize',20,'LineWidth',3,'Color',[0.85 0.33 0.1 0.5]),
    ylim([0 yaxismax])
    set(gca,'XTIck',1:3:months_tot,'XTickLabel', datestr(xbin(1:3:end-1),'yy-mmm'),'XTickLabelRotation',90)
    ymax = get(gca,'ylim');
    rectangle('Position',[1 0 14 ymax(2)],'FaceColor',[0 0 0 0.3])
end
%% plots of hours per task
yaxismax = 1.05*max(max(tasks_laborhours_by_month_by_symptom./tasks_counts_by_month_by_symptom));

figure
plot(sum(tasks_laborhours_by_month_by_symptom,2)./sum(tasks_counts_by_month_by_symptom,2),'.-','MarkerSize',28,'LineWidth',3);%(bar(sum(h.Values,2)),
set(gca,'XTIck',1:months_tot,'XTickLabel', datestr(xbin(1:end-1),'yy-mmm'),'XTickLabelRotation',90)
%set(gca,'XTIck',1:3:months_tot,'XTickLabel', datestr(xbin(1:3:end-1),'yy-mmm'),'XTickLabelRotation',90)
ylim([0 yaxismax])
ymax = get(gca,'ylim');
rectangle('Position',[1 0 14 ymax(2)],'FaceColor',[0 0 0 0.3])
title('Truck Labor Hours per Task by Month')

figure,annotation('textbox',[0.4 1 0.2 0.05],'String','Labor Hours per Task by Month',...
    'FitBoxToText','off','LineStyle','none');
for n=1:length(ybin)-1
    yaxismax = 1.1*max(tasks_laborhours_by_month_by_symptom(:,n)./tasks_counts_by_month_by_symptom(:,n));

    subplot(5,4,n),
    toplot = tasks_laborhours_by_month_by_symptom(:,n)./tasks_counts_by_month_by_symptom(:,n);
    plot(toplot,'.-','MarkerSize',20,'LineWidth',3)
    title([symptom_labels{n}])
    ylim([0 yaxismax])
    xlim([0 months_tot+1])
    set(gca,'XTIck',1:3:months_tot,'XTickLabel', datestr(xbin(1:3:end-1),'yy-mmm'),'XTickLabelRotation',90)
    %ymax = get(gca,'ylim');
    length(toplot)
    rectangle('Position',[1 0 14 yaxismax],'FaceColor',[0 0 0 0.3])
    
end
