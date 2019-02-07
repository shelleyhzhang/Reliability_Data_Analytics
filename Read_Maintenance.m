function Read_Maintenance
% Author: Shelley Zhang
% Email:  szhang@ltk.com
% Date:   May 2017
clear all
projfolder = 'P:\C4314 MBTA-X29PS01 Evaluation of Light Rail Vehicles\Type 8 Reliability Analysis\MCRS\';
projfolder = 'C:\Users\SZhang\Documents\MATLAB\T8\';
datafolder = [projfolder '\DATA_mat'];
filelist = ls([datafolder '\Mileages_2015-07-01_to_2018-12-31.mat']);
filenumb = size(filelist,1);
for n = 1%:filenumb
    load([datafolder '\' filelist(n,:)])
end
units = length(Mileage.Unit);
folder = [projfolder '\VehicleHistory\'];
filelist = ls([folder 'Vehicle History*.xls'])
files = size(filelist,1);
Maintenance = struct;
%%%%%%%%%%%%%%%%%
%Maintenance.start_date = datenum('Jul-01-2015');
Maintenance.start_date = Mileage.start_date;
%Maintenance.end_date = datenum('Aug-15-2018');
Maintenance.end_date = Mileage.end_date;
%%%%%%%%%%%%%%%%%
Maintenance.Fieldname = ...
        {'Labor Date', 'Unit', 'Meter at the time', 'Work Order', 'Total Hours',...
        'Proj #', 'Warranty YES/NO',	'Action: Description',	'Task: Description',...
        'Employee or Vendor Entered by', 'Labor Comments',	'Employee or Vendor Entered by 2nd','EMP Comments','Addtional Notes'};
% To claim memory for this Incident.Content, estimate 1500 incidents per year
%Maintenance.Content = cell(15000,12);
% Last 3 years should have 3x4=12 inspections
units_inspection_dates = cell(units,15);
Maintenance.Counter = 0;
%% Store Inspection Dates under Individual Unit
%initialize the cell string as 0
sysnum_maint = 12; %based on A,B,C,D,E,H,J,M,P,S underT,T,Y,Z.
% Letter for System under Maintenance 
sys_mainten = ['ABCDEHJMPTYZ']';%cellfun(@(c)c(1),sys_label_mainten);
DateInspect_by_Unit_o = mat2cell(zeros(units,1),ones(units,1));
MileInspect_by_Unit_o = mat2cell(zeros(units,1),ones(units,1));

DateMainten_by_Unit_Sys_o = mat2cell(repmat(0,units,sysnum_maint),ones(units,1),ones(sysnum_maint,1));
MileMainten_by_Unit_Sys_o = mat2cell(repmat(0,units,sysnum_maint),ones(units,1),ones(sysnum_maint,1));

for n = 1:files
    if n>1, clear raw txt num ,end
    %num excludes top and bottom rows compared to txt and raw
    [num, txt, raw] = xlsread([folder filelist(n,:)]);% [~,~,raw] allows nulling num and text

    %Search empty cell in the First column of txt
    rows_empty = find(cellfun(@isempty,txt(1:end-3,1)));
    %locates beginning rows of each workorder entry
    rows_entry = find(~isnan(cell2mat(raw(2:end-3,3))))+1; 
    %comment rows are supplementary to empty and entry rows
    rows_comment = setdiff([2:size(raw,1)-3],[rows_empty;rows_entry])';
    %rows that list unit number, maybe empty in older files
    rows_unit = find(strncmp(txt(:,1),'038',3));
    %inspection rows counts task code Y-8M
    rows_inspect = find(strcmp(raw(:,9),'Y-8M:TYPE 8 INSPECTION'));
    %maintenance rows of all task codes except those under Y-
    rows_mainten = setdiff(rows_entry,find(strncmp(raw(:,9),'Y-',2)));
    
    %find employee comment lines and store the content
    emp_comment = cell(length(rows_entry),1);
    rows_emp_comment = find(strncmp(txt(:,1),'EMP #',5));
    ind_in_entry = find(strncmp(raw(rows_entry+1,1),'EMP #',5));
    emp_comment(ind_in_entry) = raw(rows_emp_comment,1);
    
    %find additional notes lines and store the content
    additional_notes = cell(length(rows_entry),1);
    rows_start_additional_notes  =find(strncmp(txt(:,1),'Additional',10));
    rows_additional_notes = setdiff([2:size(raw,1)]',[rows_empty;rows_entry;rows_emp_comment;rows_unit]);
    for r=1:length(rows_start_additional_notes)
        ind_in_entry = ismember(rows_entry,...
            [max(rows_empty(rows_empty<rows_start_additional_notes(r)-2))+1:rows_start_additional_notes(r)-2]);
        tmp2 = [rows_entry;size(raw,1)-3];%in case additional notes are at the end of the file
        tmp = strcat(raw(rows_start_additional_notes(r):...
            min(tmp2(tmp2>rows_start_additional_notes(r)))-3),'\\');
        additional_notes(ind_in_entry) = {strcat(tmp{:})};
    end    
    if n>1
        Maintenance.Content = [Maintenance.Content;[raw(rows_entry,:) emp_comment additional_notes]];
    else
        Maintenance.Content = [raw(rows_entry,:) emp_comment additional_notes];
    end
    Maintenance.Counter = Maintenance.Counter + length(rows_entry);
    %unique WorkOrder on 4th column
    [~,ia,~]=unique(raw(rows_inspect,4),'stable');
    ia = rows_inspect(ia);
    %obtain the list of units under these workorders
    Inspection_Unit = str2num(cell2mat(raw(ia,2)));
    %Allocate datenum (33% less memory space than datetime) to corresponding unit
    for u = 1:length(Inspection_Unit)
        %append the new datenum under the corresponding unit
        DateInspect_by_Unit_o{ismember(Mileage.Unit,Inspection_Unit(u))}=...
            ([DateInspect_by_Unit_o{ismember(Mileage.Unit,Inspection_Unit(u))} datenum(raw{ia(u),1})]);
        MileInspect_by_Unit_o{ismember(Mileage.Unit,Inspection_Unit(u))}=...
            ([MileInspect_by_Unit_o{ismember(Mileage.Unit,Inspection_Unit(u))} raw{ia(u),3}]);
    end
    % Assign S-BA Suspension under Truck system
    Taskcode_Mainten = char(strrep(raw(rows_entry,9),'S-BA:SUSPENSION','T-S-BA:SUSPENSION'));
    % Collect maintenance system names, each file may not contain all
    if n==1
        sys_label_mainten = unique(raw(contains(txt(:,9),'-:'),9));
    else %if length(mainten_sys_label)<12
        sys_label_mainten = union(sys_label_mainten,unique(raw(contains(txt(:,9),'-:'),9)));
    end

    % Concatenate WorkOrder and System as each Workorder maybe cover several system codes
    Workorder_Sys = strcat(raw(rows_entry,4),Taskcode_Mainten(:,1));
    % Unique ID combining WorkOrder and System
    [~,ia_all,~] = unique(Workorder_Sys,'stable');
    % Rows for unique WorkOrder-System ID
    rows_unique_workorder_sys = rows_entry(ia_all);
    % System names under each unique WorkOrder-System
    Sys_list = (Taskcode_Mainten(ia_all,1));
    % Grouping by Unit number and System
    [G,ID1,ID2] = findgroups(raw(rows_unique_workorder_sys,2),cellstr(Sys_list));
    % Find Index of sorted groups
    [~,idG] = sort(G);
    % Allocate dates (sorted according to idG) into cell array, 
    % the length of each element is according to histogram counts of group IDs (G).
    Date_grouped = mat2cell(raw(rows_unique_workorder_sys(idG),1),histcounts(G,0.5:max(G)+0.5));
    Mile_grouped = mat2cell(raw(rows_unique_workorder_sys(idG),3),histcounts(G,0.5:max(G)+0.5));
    % Locate unit list ID1 in Mileage.Unit
    [~,Locunit] = ismember(str2double(ID1),Mileage.Unit);   
    % Locate system list in system under maintenance
    [~,Locsys] = ismember(ID2,sys_mainten);
%   It does not work----DateMainten_by_Unit_Sys (Locunit,Locsys)= Date_grouped; 
    for ing = 1:max(G) %for loop in groups
        %Each G corresponds to a unique combination of unit and system,
        %drop the dates in that group (ing) into the corresponding Locunit and Locsys
        DateMainten_by_Unit_Sys_o{Locunit(ing),Locsys(ing)} = ...
            [DateMainten_by_Unit_Sys_o{Locunit(ing),Locsys(ing)};Date_grouped{ing}];
        MileMainten_by_Unit_Sys_o{Locunit(ing),Locsys(ing)} = ...
            [MileMainten_by_Unit_Sys_o{Locunit(ing),Locsys(ing)};Mile_grouped{ing}];

    end
    fprintf('Completed Reading %2g:%s\n',n, filelist(n,:))
end 
%%%%%%%%%%%%%%%%%
%to capture work orders in the date range
rows = find(datenum(Maintenance.Content(:,1))>Maintenance.end_date | datenum(Maintenance.Content(:,1))<Maintenance.start_date);
Maintenance.Content(rows,:)=[];
alldates = datetime(Maintenance.Content(:,1));
% Maintenance.start_date = min(alldates);%Maintenance.start_date = Incident.start_date;
% Maintenance.end_date = max(alldates);%Maintenance.end_date = Incident.end_date;
%%check continuity
unidates = unique(alldates);
if sum(diff(unidates)>1)
    display('Missing dates between')
    missid = find(diff(unidates)>1);
    unidates([missid missid+1])
end
%%%%%%%%%%%%%%%%%
%% Grouping Maintenance dates under 12 System by 95 Unit bins
%%%%%%To update in different project%%%%%%%%%
%months_tot = months(Maintenance.start_date,Maintenance.end_date+1); 
months_tot = calmonths(between(datetime(datestr(Maintenance.start_date)),datetime(datestr(Maintenance.end_date+1))));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%end_date = addtodate(start_date,months_tot,'month');
close all
opt=4
for u = 1:units
    % clean up same inspecton within same date but under different workorders dates due to manual entry errors (+/- hours)
    % uniquetol also sorts the date, sets tolerance to 10 days or 1000 miles
    [DateInspect_by_Unit{u}, IA, IC] = uniquetol(DateInspect_by_Unit_o{u}(2:end),10/max(DateInspect_by_Unit_o{u}));
    MileInspect_by_Unit{u} = sort(MileInspect_by_Unit_o{u}((IA+1)'));
    for sysloop = 1%:sysnum_maint
        % First dummy date is 0. unique will exclude same-day workorder-system and also sort dates
         DateMainten_by_Unit_Sys{u,sysloop} = unique(datenum(DateMainten_by_Unit_Sys_o{u,sysloop}(2:end)));
         MileMainten_by_Unit_Sys{u,sysloop} = uniquetol(cell2mat(MileMainten_by_Unit_Sys_o{u,sysloop}(2:end)));
         MaintenCt_by_Unit_Sys(u,sysloop) = length(DateMainten_by_Unit_Sys{u,sysloop});
         if u==1 size(DateMainten_by_Unit_Sys), end
         switch opt
             case 1
                 % Plan A: to plot maintenance dates by Time-along X and Workorder Occurrence-along Y
                 figure(sysloop),
                 subplot(2,1,1),hold on,h = plot(DateMainten_by_Unit_Sys{u,sysloop},1:length(DateMainten_by_Unit_Sys{u,sysloop}),'.-','MarkerSize',21);
                 subplot(2,1,2),hold on,m = plot(MileMainten_by_Unit_Sys{u,sysloop},1:length(MileMainten_by_Unit_Sys{u,sysloop}),'.-','MarkerSize',21);
             case 2
                 % Plan B: to plot maintenance dates by Time-along X and differnet altitude by Unit 
                 figure(sysnum_maint+sysloop),hold on,
                 subplot(2,1,1),h2 = plot(DateMainten_by_Unit_Sys{u,sysloop},ones(1,length(DateMainten_by_Unit_Sys{u,sysloop}))*u,'.-',...
                 'MarkerSize',21,'MarkerEdgeColor','k');%,'MarkerEdgeColor',[0.2 0.5 0.2]);
                 subplot(2,1,2),m2 = plot(MileMainten_by_Unit_Sys{u,sysloop},ones(1,length(MileMainten_by_Unit_Sys{u,sysloop}))*u,'.-',...
                 'MarkerSize',21,'MarkerEdgeColor','k');%,'MarkerEdgeColor',[0.2 0.5 0.2]);
             case 3
                 %Plan C: to plot 3D
                 figure(sysnum_maint*2+sysloop),hold on,
                 subplot(2,1,1),h3 = plot3(DateMainten_by_Unit_Sys{u,sysloop},ones(1,length(DateMainten_by_Unit_Sys{u,sysloop}))*u,...
                 1:length(DateMainten_by_Unit_Sys{u,sysloop}),'.-','MarkerSize',21,'MarkerEdgeColor','k');
                 subplot(2,1,2),m3 = plot3(MileMainten_by_Unit_Sys{u,sysloop},ones(1,length(MileMainten_by_Unit_Sys{u,sysloop}))*u,...
                 1:length(MileMainten_by_Unit_Sys{u,sysloop}),'.-','MarkerSize',21,'MarkerEdgeColor','k');
             case 4
                 %disp('no plotting')
         end
        if ~isempty(DateMainten_by_Unit_Sys{u,sysloop})
            switch opt
                case 1
                    % Plan A text label
                    figure(sysloop),
                    subplot(2,1,1),text(max(DateMainten_by_Unit_Sys{u,sysloop})+2,length(DateMainten_by_Unit_Sys{u,sysloop}),...
                    num2str(Mileage.Unit(u)),'Color',h.Color,'FontSize',12)
                    subplot(2,1,2),text(max(MileMainten_by_Unit_Sys{u,sysloop})+2,length(MileMainten_by_Unit_Sys{u,sysloop}),...
                    num2str(Mileage.Unit(u)),'Color',m.Color,'FontSize',12)
                case 2
                    % Plan B text label
                    figure(sysnum_maint+sysloop),
                    subplot(2,1,1),text(mod(u,2)*(max(DateMainten_by_Unit_Sys{u,sysloop})+5)+mod(u+1,2)*(min(DateMainten_by_Unit_Sys{u,sysloop})-50),...
                    u, num2str(Mileage.Unit(u)),'Color',h2.Color)
                    subplot(2,1,2),text(mod(u,2)*(max(MileMainten_by_Unit_Sys{u,sysloop})+5)+mod(u+1,2)*(min(MileMainten_by_Unit_Sys{u,sysloop})-50),...
                    u, num2str(Mileage.Unit(u)),'Color',m2.Color)
                case 3
                    % Plan C text label
                    figure(sysnum_maint*2+sysloop),
                    subplot(2,1,1),text(max(DateMainten_by_Unit_Sys{u,sysloop})+2,u,length(DateMainten_by_Unit_Sys{u,sysloop}),...
                        num2str(Mileage.Unit(u)),'Color',h3.Color,'FontSize',12)
                    subplot(2,1,2),text(max(MileMainten_by_Unit_Sys{u,sysloop})+2,u,length(MileMainten_by_Unit_Sys{u,sysloop}),...
                        num2str(Mileage.Unit(u)),'Color',m3.Color,'FontSize',12)
                case 4
                    %disp('no plotting')
            end
        end
        % At last step, X axis label date strings with 3-month intervals
        if u == size(DateMainten_by_Unit_Sys_o,1)
            if opt<4
                subplot(2,1,1),set(gca,'XTick',datenum(datetime(datestr(Maintenance.start_date))+calmonths(0:3:months_tot)))
                subplot(2,1,1),set(gca,'XTickLabel',datestr(get(gca,'XTick')),'XTickLabelRotation',90)
                title(sys_label_mainten{sysloop})
            end
            switch opt
                case 1%Plan A
                    subplot(2,1,1),ylabel('Workorder [Count]'),xlabel('Time')
                    subplot(2,1,2),ylabel('Workorder [Count]'),xlabel('Mileage')
                case 2
                    subplot(2,1,1),ylabel('Car [Unit]'),yticklabels([])
                    subplot(2,1,2),ylabel('Car [Unit]'),yticklabels([])
                case 3
                    subplot(2,1,1),ylabel('Car [Unit]'),yticklabels([]),zlabel('Workorder [Count]')
                    subplot(2,1,2),ylabel('Car [Unit]'),yticklabels([]),zlabel('Workorder [Count]')
                case 4
                    %disp('no plotting')
%             view(2)%x as time, y as car unit
%             view(3)%rotated 3D view
%              view(0,0)%x as time, y as workorder count
            %CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10],'WellMadeVid',OptionZ)
%             OptionZ.FrameRate=10;OptionZ.Duration=11;OptionZ.Periodic=false;
%             CaptureFigVid([0,90;0,0;90,0],'WellMadeVid',OptionZ)
            end

        end
    end
end
save([folder '\..\DATA_mat\VehicleHistory_' datestr(Maintenance.start_date,'yyyy-mm-dd') '-' datestr(Maintenance.end_date,'yyyy-mm-dd') '.mat'],...
    'Maintenance','DateInspect_by_Unit','MileInspect_by_Unit','DateMainten_by_Unit_Sys','MileMainten_by_Unit_Sys','sys_mainten','sys_label_mainten')
%figure,hold on; cellfun(@plot,DateInspect_by_Unit); %only lines
%% Identify unique work order and system for Task Counts

%% Bar per Unit
close all
screensize=get(0,'ScreenSize')
figure('Position',screensize)
maxtmp = ceil(max(MaintenCt_by_Unit_Sys(:))/10)*10;
for u = 1:units %Incident
    [sorted,sorti] = sort(MaintenCt_by_Unit_Sys(u,:),'descend');
    bar(sorted,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)
    %colormap(flipud(summer(256)))
    xticks([1:size(sys_mainten,1)]);
    set(gca,'XTickLabel',sys_label_mainten(sorti),'XTickLabelRotation',30,'FontSize',14,'FontWeight','bold')
    ylabel('Workorder Tasks [Count]','FontSize',20,'FontWeight','bold')
    title(['Unit 0' num2str(Mileage.Unit(u)) ' Maintenance Task Counts by System ' datestr(Maintenance.start_date) ' to ' datestr(Maintenance.end_date)])
    xlim([0 size(sys_mainten,1)+1])
    ylim([0 maxtmp])
    saveas(gcf,[projfolder '..\figs\UnitFigs\' num2str(Mileage.Unit(u)) '\' num2str(Mileage.Unit(u)) '_MaintenSys_Bargraph'...
        datestr(Maintenance.start_date) 'to' datestr(Maintenance.end_date) '.png'])
end
%% Average Inspection gap in Days and Miles
xd = cell2mat(DateInspect_by_Unit);
xm = cell2mat(MileInspect_by_Unit);
yd = ones(size(xd,1),1)*(1:size(xd,2));
ym = yd;%ones(size(xm,1),1)*(1:size(xm,2));
% figure,subplot(2,1,1),plot(yd,xd,'o'),title('date between inspection'),xtitle('
% subplot(2,1,2),plot(ym,xm,'o'),title('mile inspect')
figure,
for u=1:units
    subplot(2,1,2),hold on,
    subplot(2,1,1),hold on,
    if ~isempty(DateInspect_by_Unit{u})
        plot(ones(length(DateInspect_by_Unit{u}),1)*u,DateInspect_by_Unit{u},'.-','MarkerSize',16,'Color',[0.2 0.5 0.2])
        subplot(2,1,2),plot(ones(length(MileInspect_by_Unit{u}),1)*u,MileInspect_by_Unit{u},'.-','MarkerSize',16,'Color',[0.2 0.5 0.2])
        d{u} = diff(DateInspect_by_Unit{u});
        m{u} = diff(MileInspect_by_Unit{u});
    end
end
subplot(2,1,1),
ylabel('Inspection Dates [Days]','FontSize',18,'FontWeight','bold')
ylim([min(xd)-1 max(xd)+1]);
ry = get(gca,'YTick');yticklabels(ry-min(xd)+1)
xlabel('Type 8 Car Unit','FontSize',18,'FontWeight','bold')
set(gca,'XTick',[1:units]);xlim([1 units])
set(gca,'XTickLabel',Mileage.Unit,'XTickLabelRotation',90,'FontWeight','bold')

subplot(2,1,2),
ylabel('Inspection Mileage [Miles]','FontSize',18,'FontWeight','bold')
xlabel('Type 8 Car Unit','FontSize',18,'FontWeight','bold')
set(gca,'XTick',[1:units]);xlim([1 units])
set(gca,'XTickLabel',Mileage.Unit,'XTickLabelRotation',90,'FontWeight','bold')

%All Inspection Gaps in days 
DateInspect_gaps = cell2mat(d);
MileInspect_gaps = cell2mat(m);
%DateInspect_gaps(DateInspect_gaps>120)=[];
figure,histogram(DateInspect_gaps,min(100,range(DateInspect_gaps))),set(gca,'YScale','log')
figure,histogram(MileInspect_gaps,min(100,range(MileInspect_gaps)))
[min(nonzeros(DateInspect_gaps)) max(nonzeros(DateInspect_gaps)) mean(nonzeros(DateInspect_gaps)) std(nonzeros(DateInspect_gaps))]
[min(nonzeros(MileInspect_gaps)) max(nonzeros(MileInspect_gaps)) mean(nonzeros(MileInspect_gaps)) std(nonzeros(MileInspect_gaps))]
% [mean(nonzeros(m)) std(nonzeros(m))]
%% Failure in Service dates grouping by 95 units to count Mean Time to Failure from last inspection (up to 90 days)
%% Failure in Service dates grouping by 95 units and 11 Systems, to count MTF for each system
%close all
%load([folder '\..\DATA_mat\VehicleHistory' datestr(Maintenance.start_date) '-' datestr(Maintenance.end_date) '.mat'])
%Incident.Content has column 1,3,5 as unit, date, symptom code
Incidents_Counted = Incident.Content(Incident.events_counted,[1 3 5]);
%Incidents_Counted = Incident.Content(:,[1 3 5]);
Incidents_by_unit = discretize(str2double(Incidents_Counted(:,1)),[Mileage.Unit-0.5 max(Mileage.Unit)+0.5]);
Incidents_by_sys = findgroups(double(incident_by_system(Incident.events_counted)));
sys_incident = unique(incident_by_system(Incident.events_counted));

sys_compare = 'ABCDEHJPT';%9 systems to compare between incidents and maintenance
sys_compare_incident = ismember(sys_incident,sys_compare);
sys_compare_mainten = ismember(sys_mainten,sys_compare)';
%initialize to NaN
MeanDays_from_inspection = NaN(units,1);
%MeanMiles_from_inspection = NaN(units,1);
MeanDays_from_mainten = NaN(units,length(sys_compare));
%MeanMiles_from_mainten = NaN(units,length(sys_compare));
for u = 1:units
    if ~isempty(find(Incidents_by_unit==u))&& ~isempty(DateInspect_by_Unit{u})
        %sort the incident dates with matching unit as u
        DateIncident_by_unit = sort(datenum(Incidents_Counted(find(Incidents_by_unit==u),2)));
      
        %+1 to count starting from the day following maintenance since
        %incidents can happen on the same day to induce maintenance rather than after maintenance
        DateIncident_between_inspection = discretize(DateIncident_by_unit',[DateInspect_by_Unit{u}+1 max(DateInspect_by_Unit{u})+90]);
        %DateIncident_between_inspection = discretize(DateIncident_by_unit',[DateInspect_by_Unit{u} max(DateInspect_by_Unit{u})+90],'IncludeEdge','right');
        DateIncident_between_inspection(isnan(DateIncident_between_inspection))=-1;
        %find the earliest date after inspection
        [C,ia,ic] = unique(DateIncident_between_inspection);
        %because dates_between_inspection contains -1 (converted from NaN),
        DaysIncident_from_inspection{u} = round([DateIncident_by_unit(ia(2:end))'-DateInspect_by_Unit{u}(C(2:end))]);
        MeanDays_from_inspection(u) = round(mean(DaysIncident_from_inspection{u}));
    end
    for sysloop = 6%1:length(sys_compare)
            si = find(sys_incident==sys_compare(sysloop));
            sm = find(sys_mainten ==sys_compare(sysloop));
            %sort the incident dates with matching unit as u and system as s
            usi = intersect(find(Incidents_by_unit==u),find(Incidents_by_sys==si));
            if ~isempty(usi)
                DateIncident_by_unit_sys = sort(datenum(Incidents_Counted(usi,2)));
                if ~isempty(DateIncident_by_unit_sys)&& ~isempty(DateMainten_by_Unit_Sys{u,sm})
                    DateIncident_between_mainten = discretize(DateIncident_by_unit_sys,[DateMainten_by_Unit_Sys{u,sm}+1 ;Maintenance.end_date+1]);
                    %convert NaN to some number, -1, which is distinctive from all other possibilities
                    DateIncident_between_mainten(isnan(DateIncident_between_mainten))=-1;
                    %find the earliest date after maintenance
                    [C2,ia2,ic2] = unique(DateIncident_between_mainten);
                    DaysIncident_from_mainten = floor([DateIncident_by_unit_sys(ia2(2:end))-DateMainten_by_Unit_Sys{u,sm}(C2(2:end))]);
                    %[u; s; nonzeros(Days_from_mainten)]
                    MeanDays_from_mainten(u,sysloop) = round(mean(nonzeros(DaysIncident_from_mainten)));
                    if u==units
                        my_bar_avgstd(MeanDays_from_mainten(:,sysloop),Mileage.Unit)
                        title(['Mean Time to Failure from Maintenance under ' sys_label_mainten{sm}(4:end)])    
                    end
                else
                    if u==units
                        my_bar_avgstd(MeanDays_from_mainten(:,sysloop),Mileage.Unit)
                        title(['Mean Time to Failure from Maintenance under ' sys_label_mainten{sm}(4:end)])    
                    end
                end
            else
                display(['Unit ' num2str(UNIT(u)) '-sys' sys_compare(sysloop) ' has empty intersection'])

                if u==units
                        my_bar_avgstd(MeanDays_from_mainten(:,sysloop),Mileage.Unit)
                        title(['Mean Time to Failure from Maintenance under ' sys_label_mainten{sm}(4:end)])    
                end
            end
       end
 end
my_bar_avgstd(MeanDays_from_inspection,Mileage.Unit,'hor')
%[sorted,sortedi]=sort(MeanDays_from_inspection);
title('Mean Time to Failure from Last Inspection (up to 90 days)')
%% Work Order Counts by car
%Car Unit in Incident.Content(:,1)
clear dates_inc validcars

noY = find(~startsWith(Maintenance.Content(:,9),'Y-'));
[workorder_noY,ia,~] = unique(Maintenance.Content(noY,4));
%capture work orders without Inspection Y
dates_inc = datestr(Maintenance.Content(noY(ia),1));
dates_inc = datetime(dates_inc);

workorder_in_quarter = my_discretize(datetime(datestr(Maintenance.start_date)),dates_inc,'quarter');    
quarters_inc = max(workorder_in_quarter);
workorder_in_month = my_discretize(datetime(datestr(Maintenance.start_date)),dates_inc,'month');    
months_inc = max(workorder_in_month);
monthfiller = zeros(quarters_inc*3-months_inc,1)';


cars_with_workorder_noY = str2double(unique(Maintenance.Content(noY(ia),2)));
validcars = ismember(Mileage.Unit,cars_with_workorder_noY);
workorder_by_car = str2double(Maintenance.Content(noY(ia),2));
%initialize
workorder_noY_by_car_by_quarter = zeros(units,quarters_inc);
workorder_noY_by_car_by_quarter(validcars,:) = histcounts2(workorder_by_car,workorder_in_quarter,...
    [cars_with_workorder_noY-0.5;(cars_with_workorder_noY(end)+0.5)], [0.5:1:(quarters_inc+0.5)] );
    
%%%%%%%%%%%%%Out of Service Cars%%%%%%%%%%%%%%%%
cars_oos = [3803 3807 3808 3832 3836 3854 3873 3879 3891 3893]; %3879 not even show up in Inspections Due
% validcars(ismember(Mileage.Unit,cars_oos))=0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mean Miles between Work Orders per Quarter, excluding work orders solely for Y- inspection tasks.
firstday_in_quarter = datenum(datetime(datestr(Maintenance.start_date))+calquarters(0:quarters_inc));
miles_by_car_by_quarter = zeros(units,quarters_inc);
%since Mileage was since 2010, here we just capture the last 

miles_by_car_by_quarter = squeeze(sum(reshape([Mileage.Vals(:,size(Mileage.Vals,2)-months_inc+1:end) repmat(monthfiller,units,1)],units,3,quarters_inc),2));

tmp = workorder_noY_by_car_by_quarter;
tmp(tmp==0)=1; 
MMBWO_car_quarter = miles_by_car_by_quarter./tmp;
[B,I]=sort(MMBWO_car_quarter(:,end),'descend');
%display the latest quarter in rank of decreasing MMBWO
[Mileage.Unit(I)' miles_by_car_by_quarter(I,end) workorder_noY_by_car_by_quarter(I,end) B]

pause
[B,I]=sort(workorder_noY_by_car_by_quarter(:,end));
cars_oos = [3807        3854        3873        3879        3891        3808        3832        3803        3893];
tmp = [Mileage.Unit(I)'  B];
tmp (ismember(I,cars_oos-3799),:)=[]
%% Fleet Labor Hours and Task Counts by Quarter by System
tasks_by_system = char(Maintenance.Content(noY,9));
tasks_by_system = tasks_by_system(:,1);
% tasks_by_system = mat2cell(tasks_by_system(:,1),ones(size(tasks_by_system,1),1),1);
%to map char numbers 65, 66,...90 to continuous number 1 to 12
tmp = unique(double(tasks_by_system));
sysnum(tmp)=[1:length(tmp)];
tasks_by_sysnum = sysnum(double(tasks_by_system))';

tasks_by_quarter = ((categorical(discretize(datetime(Maintenance.Content(noY,1)),'quarter','categorical'))));
[Group,IDs,IDq] = findgroups(tasks_by_sysnum,double(tasks_by_quarter));
tasks_laborhours_by_quarter_by_system = zeros(max(IDq),max(IDs));
newsub = sub2ind([max(IDq) max(IDs)],IDq,IDs);
tmp = splitapply(@sum,cell2mat(Maintenance.Content(noY,5)),Group);%returns a row
tasks_laborhours_by_quarter_by_system(newsub) = tmp;

%locate unique combination of work order and task system
[C,ia,ic] = unique(join(Maintenance.Content(noY,[4 9])));
[C2,ia,ic] = unique([Maintenance.Content(noY,4),mat2cell(tasks_by_system,ones(size(tasks_by_system)),1)]);
[C3,ia,ic] = unique(join([Maintenance.Content(noY,4),mat2cell(tasks_by_system,ones(size(tasks_by_system)),1)]));

%% Mean Time to Failure from Last Maintenance by System
% systems = length(system);
% for n=1:size(Incidents_Counted,1)
%     incidents_counted_by_system(n) = Incidents_Counted{n,3}(1);
% end
% system = unique(incidents_counted_by_system); %'ABCDEHJMPTVZ'
% systems = length(system);
% for s = 1:system
    
    
    
%%
%     sz = size(raw); %should be 4n+3 by 8
%     blockrow = 4;
%     blockcol = sz(2);
%     Incident.Counter(n+1) = uint16((sz(1)-3)/blockrow);
%     pickcell = [1:14 17 25]'+1*blockcol; %the 16 cells with content
%     %tmp = reshape(reshape(raw(2:end-2,:)',[8 blockrow (sz(1)-3)/blockrow]),[8*blockrow (sz(1)-3)/blockrow])';
%     picklocs = uint16(linspaceNDim(pickcell,pickcell+(Incident.Counter(n+1)-1)*blockrow*blockcol,Incident.Counter(n+1)));
%     raw = raw';
%     Incident.Content((1+sum(Incident.Counter(1:n))):sum(Incident.Counter(1:n+1)),:) = raw((picklocs))';
%     fprintf('Read %2g:%s\n',n, filelist(n,:))


