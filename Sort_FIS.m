function Sort_FIS
% Author: Shelley Zhang
% Email:  szhang@ltk.com
% Date:   Oct 2017
% Phase II
clear all
folder = 'P:\C4314 MBTA-X29PS01 Evaluation of Light Rail Vehicles\Type 8 Reliability Analysis\MCRS\DATA_mat';
folder = 'E:\DATA_mat';
folder = 'C:\Users\SZhang\Documents\MATLAB\T8\DATA_mat';
filelist = ls([folder '\*_2015-07-01_to_2018-12-31.mat']);
filenumb = size(filelist,1);
for n = 1:filenumb
    load([folder '\' filelist(n,:)])
end
units = length(Mileage.Unit);
system = unique(incident_by_system);
systems = length(system);
start_date = datetime(datestr(Incident.start_date));%datetime('01-Apr-2010');
start_year = (squeeze(start_date.Year))+1; %fiscal year 2011 is 2010-Jul to 2011-Jun
%% Histcounts Incidents by System in 5th column, contains A,B,C,D,E,H,J,M,P,T,Z
[~,incident_by_car]=ismember(str2num(cell2mat(Incident.Content(:,1))),Mileage.Unit);
dates = datestr(Incident.Content(:,3));
dates = datetime(dates);
incident_in_month = discretize(dates,'month');%last 2 years contains 61-84 months
months = max(incident_in_month);
incident_in_quarter = discretize(dates,'quarter'); 
quarters = max(incident_in_quarter);
incident_in_year = ceil(incident_in_quarter/4);
year = unique(incident_in_year)+start_year-1;
years = max(incident_in_year);
system_label = {'A (Air Supply)','B (Brake)','C (Vehicle Body & Frame)','D (Doors)','E (Electrical and Lighting)','H (HVAC)',...
                'J (Coupler and Draft Gear)','P (Engine)', 'T (Truck)','V (Automatic Fare Collection)','Z (Other)'};
system_in_double = double(system(:));
incident_counts_by_system = histcounts(double(incident_by_system),[system_in_double-0.5; max(system_in_double)+0.5]);
incident_count_by_system_year = histcounts2(double(incident_by_system), incident_in_year',...
                                 [system_in_double-0.5; max(system_in_double)+0.5],0.5:(years+0.5));
incident_count_by_system_quarter = histcounts2(double(incident_by_system), incident_in_quarter',...
                                 [system_in_double-0.5; max(system_in_double)+0.5],0.5:(quarters+0.5));
                            
filler  = zeros(years*4-quarters, 1);
%miles_by_year = sum(reshape([filler sum(MileageMile,1)],4,years),1);   
%%%%%%%%%%%%%%%%%%%%
incident_count_by_system_quarter'
pause
%%%%%%%%%%%%%%%%%%%%
%% Incident by Subsystem (Symptom )by Car over 7 years 
close all
symptom_counter = 0;
for n = 1:systems
    
    ind = findstr(incident_by_system,system(n));
    symptom_under_system_n = Incident.Content(ind,5);
    symptoms_under_system_n = length(unique(symptom_under_system_n));
    
    date_under_system_n = Incident.Content(ind,3);
%     switch n
%         case 1
%             symptom_under_system_n(contains(symptom_under_system_n,{'A01','A02','A09','A10'}))={'A01-02-09-10 (LOW/NO AIR)'};
%             symptom_under_system_n(contains(symptom_under_system_n,{'A31','A41','A42'}))={'A31-41-42 (AIR SPRING/BAG)'};
%         case 8
%             symptom_under_system_n(contains(symptom_under_system_n,{'P01','P02','P11','P27'}))={'P01-02-11-27 (LOST POWER)'};
%             symptom_under_system_n(contains(symptom_under_system_n,{'P51','P52'}))={'P51-52 (P/D FAULT INDICATION A&B)'};
%     end
    clear counts* cats* ind2
    symptomcats_under_system_n = categorical(symptom_under_system_n);
    [counts,cats] = histcounts(symptomcats_under_system_n);
    system_label{n}=char(cats(1));
    symptom_label{n} = cats;
    symtpom_label_all(symptom_counter+1 : symptom_counter+symptoms_under_system_n)=cats;   
    %pie
    figure,pie(symptomcats_under_system_n),colormap(flipud(summer(256)))
    close;
    
    %bar horizontal ranked by counts large to small    
    [counts2,ind2] = sort(counts);
    figure,barh(counts2,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)
    yticks(1:length(cats)), yticklabels(cats(ind2))
    ylim([0 length(cats)+1])
    xlabel('Occurrences 04/2010 - 03/2017 [Count]')
    close;

    %by car by symptoms
    %grp2idx only exists within Stats toolbox
    %symptomind_under_system_n = grp2idx(symptomcats_under_system_n);
    symptomind_under_system_n = double(symptomcats_under_system_n);
    
    h=histogram2(symptomind_under_system_n,incident_by_car(ind),[0.5:symptoms_under_system_n+0.5],[0.5:units+0.5]);
    symptom_dist_by_car(symptom_counter+1 : symptom_counter+symptoms_under_system_n,:) = h.Values;
    incident_counts_by_car_under_system_n (n,:) = sum(h.Values,1); 
    symptom_counter = symptom_counter + symptoms_under_system_n;
 
    %histogram of Incidents under System over Year
    figure, lastyears = min(years, 7 )
    symptom = unique(symptom_under_system_n);
    symptoms = length(symptom);
    sympnum = mat2cell(num2str([1:symptoms]'),ones(symptoms,1));
    incident_sympnum = double(categorical(symptom_under_system_n,symptom,sympnum));
    incident_quarter = my_discretize(start_date,datetime(date_under_system_n),'quarter');
    incident_year = ceil(incident_quarter/4);
    h = histcounts2(incident_sympnum,incident_year,[0.5:1:symptoms+0.5],[years-lastyears+0.5:years+0.5]);
    lh = barh(h(ind2,end-lastyears+1:end),'hist');
    title([system_label{n} ' Symptom by Year'])
    xlabel(['Occurrences [Count]'])
    neworder = [years:-1:years-lastyears+1]';
    legend(lh(end:-1:1),cellstr([repmat('Year ',[lastyears 1]) num2str(year(neworder))]),'location','bestoutside');
    yticks(1:symptoms); ylim([0 symptoms+1])
    set(gca,'YTickLabel',(symptom(ind2)))
    close
    %figure, pie(sum(h(:,end-2:end),2)),title(system_label{n})%,legend(symptom)
    %{
    lastyears=3;
    firstind = min(find(incident_in_year>years-lastyears));
    figure, 
    subplot(2,1,1),histogram(categorical(Incident.Content(ind(ind>=firstind),5))),title(system_label{n})
    subplot(2,1,2),pie(categorical(Incident.Content(ind(ind>=firstind),5))),title(system_label{n})
    %}
end
figure,imagesc(incident_counts_by_car_under_system_n)
colormap(gray)
%% Inc Counts Pie for System PerCar
close all
screensize=get(0,'ScreenSize')
figure('Position',screensize)
maxtmp = ceil(max(incident_counts_by_car_under_system_n(:))/10)*10;
for u = 1:units %Incident
    [sorted,sorti] = sort(incident_counts_by_car_under_system_n(:,u),'descend');
    bar(sorted,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)
    %colormap(flipud(summer(256)))
    xticks([1:systems]);
    set(gca,'XTickLabel',system_label(sorti),'XTickLabelRotation',30,'FontSize',14,'FontWeight','bold')
    ylabel('Incidents [Count]','FontSize',20,'FontWeight','bold')
    title(['Unit 0' num2str(Mileage.Unit(u)) ' Incident Counts by System ' datestr(Mileage.start_date) ' to ' datestr(Mileage.end_date)])
    xlim([0 systems+1])
    ylim([0 maxtmp])
    saveas(gcf,[folder '\UnitFigs\' num2str(Mileage.Unit(u)) '_IncSys_Bargraph'...
        datestr(Mileage.start_date) 'to' datestr(Mileage.end_date) '.png'])
end
%% Brake System Last two years for Leo and Nate's investigation
% T8HPCUREL is the project code for HPCU Modification
% for n = 2
%     startmonth = 1;
%     in_system_month = strcmp(cellstr(incident_by_system'),system(n))& ge(incident_in_month,startmonth);
%     figure,counts = histcounts(incident_in_month(in_system_month),startmonth-0.5:months+0.5);
%     plot(counts)
%     xticks(startmonth:months)
% 	%create a series of monthly datetime starting from 2010, 4th month for next 84 months forward
%     firstdate = min(dates);
%     category_in_month = discretize(datetime(firstdate.Year,firstdate.Month+[0:months-1],firstdate.Day),'month','categorical')
%     set(gca,'XTickLabel',char(category_in_month(startmonth:end)),'fontsize',8,'XTickLabelRotation',90)
%     %set(gca,'XTickLabelRotation',90)
%     ylabel('Occurrence [Count]','fontsize',12)
%     title('Brake System Incident By Month')
%     axis tight
%     hold on, plot(61:84,counts(61:end),'r')
% end
%% Incident Counts By System by Year
% not too helpful
% figure,imagesc((incident_count_by_system_year)./repmat(miles_by_year,[11 1]))
% yticks(1:systems); yticklabels(system_label)
% xticks(1:years); xticklabels([repmat('Year',[years 1]) num2str([1:years]')])
% title('Incident Counts by System by Year (normalized by Fleet Miles by Year)')
% colormap(hot)
yearind = [years-min(length(years)-1,6):years];
[sorted,si] = sort(sum(incident_count_by_system_year(:,yearind),2),'descend');
%figure,bar(incident_count_by_system_year(si,yearind)*1e6./repmat(miles_by_year(yearind),[systems 1]),'grouped')
figure,bar(incident_count_by_system_year(si,yearind),'grouped')
title('Incident Counts by System (normalized by Fleet Mileages in Million)')

%yticks(1:years); yticklabels([repmat('Year',[lastyears 1]) num2str(year(end-lastyears+1:end)')])
%set(gca,'xaxisLocation','top')
xticks(1:systems); xticklabels(system_label(si))
set(gca,'XTickLabelRotation',30)
legend(cellstr([repmat('Year',[length(yearind) 1]) num2str(year(yearind))]))
ylabel('Counts per Million Miles')
% figure,histogram2(double(incident_by_system), incident_in_year,...
%                                  [system_in_double-0.5; max(system_in_double)+0.5],0.5:(years+0.5))
%Incident Counts Last 3 years
yearind = [years-2:years];
[sorted,si] = sort(sum(incident_count_by_system_year(:,yearind),2));
figure,barh(sum(incident_count_by_system_year(si,yearind),2),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)%./sum(miles_by_year(yearind)
yticks(1:systems); yticklabels(system_label(si))
%set(gca,'XTickLabelRotation',90)
xlabel('Occurrence [Count]')
title('Incidents from 04/01/2014 to 03/31/2017')
text(sum(incident_count_by_system_year(si,yearind),2),1:systems, num2str(sum(incident_count_by_system_year(si,yearind),2)))
%set(gca,'TickLength',[0 0 ])
%% Per Station 9th Column
%{
%tmp = (strfind(Incident.Content(:,9),'STA*'))
incident_location = erase(Incident.Content(Incident.events_counted,9),'STATION ');
%incident_location = erase(Incident.Content(:,9),'STATION ');
incident_location = erase(incident_location,'STATIO ');
incident_location = erase(incident_location,'STATI ');
incident_location = erase(incident_location,'(B) ');
incident_location = erase(incident_location,'(B ');
incident_location = erase(incident_location,'(C) ');
incident_location = erase(incident_location,'(C ');
incident_location = erase(incident_location,'(D) ');
incident_location = erase(incident_location,'(D ');
incident_location = erase(incident_location,'(E) ');
incident_location = erase(incident_location,'(E ');
incident_location = erase(incident_location,'( ');
incident_location = erase(incident_location,'APPROACHING ');
incident_location = erase(incident_location,'APPROACH ');
incident_location = erase(incident_location,'APPR. ');
incident_location = erase(incident_location,'APPR.');
incident_location = strrep(incident_location,' T ',' ');
incident_location = strrep(incident_location,'CTR ','CENTER ');
incident_location = strrep(incident_location,'CTR. ','CENTER ');
incident_location = strrep(incident_location,'CENT ','CENTER ');
incident_location = strrep(incident_location,'GOV.','GOVERNMENT ');
incident_location = strrep(incident_location,'C ','CENTER ');
incident_location = strrep(incident_location,'KEN ','KENMORE ');
incident_location = strrep(incident_location,'RESERVIOR','RESERVOIR');
incident_location = strrep(incident_location,'EASTBOUND','INBOUND');
incident_location = strrep(incident_location,'WESTBOUND','OUTBOUND');
a=strfind(incident_location,'VILL');
for n = find(~cellfun('isempty',a))'
    tmp = strfind(incident_location(n),' ');
    incident_location{n} = ['GREEN LINE   BROOKLINE VILLAGE' incident_location{n}(tmp{1}(end):end)];
end
[counts, cats] = histcounts(categorical(incident_location));
[counts_sorted,index_sorted ] = sort(counts);
%bar graph the highest occured 20 locations
toplot = 30;
figure,barh(counts_sorted(end-toplot+1:end),'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)
ylim([0 toplot+1])
yticks(1:toplot)
yticklabels(cellfun(@(x) x(14:end), cats(index_sorted(end-toplot+1:end)), 'un', 0))
xlabel('Occurrences over 7 Years [Count]')
%% At Reservior Inbound
station_target = 'COPLEY INBOUND';
 station_target = 'RESERVOIR INBOUND';
 station_target = 'NORTH OUTBOUND';
%station_target = 'KENMORE INBOUND';
a=strfind(incident_location,station_target);
h = histcounts(categorical(cellstr(incident_by_system(find(~cellfun('isempty',a)))')),cellstr(system'));
[h2,i2] = sort(h);
figure,barh(h2,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)
yticklabels(cellstr(system_label(i2)'))
%ylim([0 systems])
xlabel('Occurrences over 7 Years [Count]')
title(['Incident by System at ' station_target])
%% To identify Type 7 vehicles (unit 36xx or 37xx) that were coupled at time of failure
% comment in 16th column
startind = regexp(Incident.Content(:,16),['3[67]\d\d\[s'],'once'); %hopefully there is only one Type7 unit in each incident
incident_type7_index = find(~cellfun(@isempty,startind));
%strindex_type7 = cell2mat(cellfun(@(c) [c(1):c(1)+3],startind(incident_type7_index),'UniformOutput',false));
%to find out the longest string at which index, the real index is incident_type_index(max_index)
%[max_size, max_index] = max(cellfun('length', Incident.Content(incident_type7_index,16)));
for n = 1:length(incident_type7_index)
    tmp = startind{incident_type7_index(n)}; %string index on the corresponding incident index
    units_type7_coupled(n) = str2num(Incident.Content{incident_type7_index(n),16}(tmp:tmp+3));
    incident_coupled_system(n) = Incident.Content{n,5}(1);
end
type7_id = unique(units_type7_coupled);
%Type 7 has units 3600-3699 and 3700-3719
type7_id = intersect([3600:3719],type7_id);
[counts] = histcounts(units_type7_coupled,[type7_id-0.5 max(type7_id)+0.5]);
[sorted,ind] = sort(counts);
figure,barh(sorted,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)
title('Type 7 Vehicles Coupled with Type 8 occurring Incidents')
yticks(1:length(type7_id))
yticklabels(type7_id(ind))
set(gca,'XTickLabelRotation',90)
xlabel('Incident Occurrences over 7 Years [Count]')

[counts,cats] = histcounts(categorical(cellstr(incident_coupled_system')));
[sorted,ind] = sort(counts,'descend');
figure,bar(sorted,'FaceColor',[0 .5 .5],'EdgeColor',[0 .9 .9],'LineWidth',1.5)
title('Systems prone to Failure when coupled with Type 7 Vehicles')
xticks(1:length(cats))
xticklabels(system_label(ind))
set(gca,'XTickLabelRotation',45)
ylabel('Occurrences over 7 Years [Count]')
%carunits_type7_coupled = cellfun(@(a,b)a(b),Incident.Content(incident_type7_index,16),strindex_type7,'uniform',0);

%% To identify Systems that are prone to failure when coupled with Type 7
%% Grouping Incident dates under 12 System by 95 Unit bins
%%%%%%To update in different project%%%%%%%%%
start_date = datenum('Apr-1-2014'); months = 36; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end_date = addtodate(start_date,months,'month');
close all
opt=4
for u = 1:units
    % clean up same inspecton within same date but under different workorders dates due to manual entry errors (+/- hours)
    % uniquetol also sorts the date, sets tolerance to 10 days or 1000 miles
    [DateI_by_Unit{u}, IA, IC] = uniquetol(DateInspect_by_Unit_o{u}(2:end),10/max(DateInspect_by_Unit_o{u}));
    MileInspect_by_Unit{u} = sort(MileInspect_by_Unit_o{u}((IA+1)'));
    for sysloop = 1:sysnum_maint
        % First dummy date is 0. unique will exclude same-day workorder-system and also sort dates
         DateMainten_by_Unit_Sys{u,sysloop} = unique(datenum(DateMainten_by_Unit_Sys_o{u,sysloop}(2:end)));
         MileMainten_by_Unit_Sys{u,sysloop} = uniquetol(cell2mat(MileMainten_by_Unit_Sys_o{u,sysloop}(2:end)));
         if u==1 size(DateMainten_by_Unit_Sys), end
         switch opt
             case 1
                 % Plan A: to plot maintenance dates by Time-along X and Workorder Occurrence-along Y
                 figure(sysloop),hold on,
                 subplot(2,1,1),h = plot(DateMainten_by_Unit_Sys{u,sysloop},1:length(DateMainten_by_Unit_Sys{u,sysloop}),'.-','MarkerSize',21);
                 subplot(2,1,2),m = plot(MileMainten_by_Unit_Sys{u,sysloop},1:length(MileMainten_by_Unit_Sys{u,sysloop}),'.-','MarkerSize',21);
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
                subplot(2,1,1),set(gca,'XTick',datenum(datetime(datestr(start_date))+calmonths(0:3:months)))
                subplot(2,1,1),set(gca,'XTickLabel',datestr(get(gca,'XTick')),'XTickLabelRotation',90)
                title(sys_label_mainten{sysloop})
            end
            switch opt
                case 1%Plan A
                    subplot(2,1,1),ylabel('Workorder [Count]')
                    subplot(2,1,2),ylabel('Workorder [Count]')
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
%}