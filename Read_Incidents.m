function Read_Incidents
% Author: Shelley Zhang
% Email:  szhang@ltk.com
% Date:   April 2017
%%
clear all
folder = 'P:\C4314 MBTA-X29PS01 Evaluation of Light Rail Vehicles\Type 8 Reliability Analysis\MCRS\MCRS_IncidentReport\';
folder = 'E:\MCRS_IncidentReport\';
folder = 'C:\Users\SZhang\Documents\MATLAB\T8\IMR\';
filelist = ls([folder '*.xls*']);
files = size(filelist,1);
Incident = struct;
%%%%%%%%%%%%%%%%%
Incident.start_date = datenum('Jul-01-2015');
Incident.end_date = datenum('Dec-31-2018');
%%%%%%%%%%%%%%%%%
Incident.Fieldname = ...
        {'CarUnit','IncID','Date','Dispatcher','FIS Symptom','Passenger Load','Service Action','Trips Lost Num',...
         'Location Info','Operator','Vehicle Action','Special Attention','Service Impact','Trip Lost Service Interruption'...
         'System Affected',...
         'Comment'};
% To claim memory for this Incident.Content, estimate 1500 incidents per year
Incident.Content = cell(15000,16);
Counter = 0; %to count incident entries in each file

for n = 1:files
    [~,~,raw] = xlsread([folder filelist(n,:)]); %can also readtable
    %raw contains both numbers and text, below 1st row, arranged in 4n by 8 blocks for each incident
    sz = size(raw); %should be 4n+3 by 8
    blockrow = 4;
    blockcol = sz(2);
    Counter(n+1) = uint16((sz(1)-3)/blockrow);
    pickcell = [1:14 17 25]'+1*blockcol; %the 16 cells with content
    %tmp = reshape(reshape(raw(2:end-2,:)',[8 blockrow (sz(1)-3)/blockrow]),[8*blockrow (sz(1)-3)/blockrow])';
    picklocs = uint16(linspaceNDim(pickcell,pickcell+(Counter(n+1)-1)*blockrow*blockcol,Counter(n+1)));
    raw = raw';
    Incident.Content((1+sum(Counter(1:n))):sum(Counter(1:n+1)),:) = raw((picklocs))';
    fprintf('Read %2g:%s\n',n, filelist(n,:))
end
%to save space when storing content
Incident.Content(sum(Counter)+1:end,:)=[];
%excludes non-equippment related incidents
Incident.Content(contains(Incident.Content(:,5),{'Z98','Z02','Z20','C01','C02','C04','C42','V07','V10','V12','V16','V17','V18'}),:)=[];
Incident.events_intrpt = strcmp(Incident.Content(:,14),'Y');%count "Service Interruption"= "Y" incidents
%to capture incidents in the date range
rows = find(datenum(Incident.Content(:,3))>Incident.end_date | datenum(Incident.Content(:,3))<Incident.start_date);
Incident.Content(rows,:)=[];
Incident.events_intrpt(rows)=[];
for n=1:size(Incident.Content,1)
    incident_by_system(n) = Incident.Content{n,5}(1);
end
towrite = [folder '\..\DATA_mat\Incidents_' datestr(Incident.start_date,'yyyy-mm-dd') '_to_' datestr(Incident.end_date,'yyyy-mm-dd')];
save([towrite '.mat'],'Incident','incident_by_system')
% xlswrite([towrite '.xlsx'],Incident.Fieldname,1,'A1')
% xlswrite([towrite '.xlsx'],Incident.Content,1,'A2')
