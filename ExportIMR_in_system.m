function ExportIMR_in_system
clear all
folder = 'C:\Users\SZhang\Documents\MATLAB\T8\DATA_mat\';
load([folder  'Incidents_01-Oct-2010_30-Oct-2017.mat'])
%for loop is faster than cellfun to extract first letter from symptom code
for n=1:size(Incident.Content,1)
    incident_by_system(n) = Incident.Content{n,5}(1);
    if Incident.Content{n,5}(2:3) == ' ('
        system_names = Incident.Content(n,5);
    end
end
system_label = unique(system_names);
%%
desired_start_date = max(datetime(2014,10,1),datetime(datestr(Incident.start_date)));
desired_end_date = min(datetime(2017,10,31),datetime(datestr(Incident.end_date)));
months_tot = calmonths(between(desired_start_date,desired_end_date+1));
%Incident has 3rd column storing dates
tf=isbetween(Incident.Content(:,3),desired_start_date,desired_end_date);
tfi = uint32(find(tf==1)); 
%Incident has 5th column storing symptom codes
event_symptom_code=char(Incident.Content(tf,5));

systems = unique(event_symptom_code(:,1));
opt = input(['System to request: D/[' systems' ']:'],'s');
if isempty(opt)
    opt = 'D';
end
if length(tfi)< 2^16
    rows1 = uint16(strfind(event_symptom_code(:,1)',opt));
    %rows where only letter space ( show up, i.e. high level system codes
    rows2 = uint16(strfind(reshape(event_symptom_code(:,2:3)',[1 length(event_symptom_code)*2]),' (')+1)/2;
else
    rows1 = uint32(strfind(event_symptom_code(:,1)',opt));
    %rows where only letter space ( show up, i.e. high level system codes
    rows2 = uint32(strfind(reshape(event_symptom_code(:,2:3)',[1 length(event_symptom_code)*2]),' (')+1)/2;
end
system_label = unique(event_symptom_code(rows2,:),'rows');
tmp = char(system_label);
system_opt = replace(strtrim(system_label(strfind(tmp(:,1)',opt),:)),' ','_');
filename = [system_opt '_Incidence_' datestr(desired_start_date,'yyyy-mm-dd') '_to_' datestr(desired_end_date,'yyyy-mm-dd') '.xlsx']
inputcon = Incident.Content(tfi(rows1),:);

xlswrite([folder filename],Incident.Fieldname)
xlswrite([folder filename],inputcon,1,'A2')
%%
incident_symptom_code=char(Incident.Content(tf,5));
systems = unique(incident_symptom_code(:,1));
opt = input(['System to request: D/[' systems' ']:']);
if isempty(opt)
    opt = 'D';
end
rows = strfind(incident_symptom_code(:,1)',opt);
%rows where only letter space ( show up, i.e. high level system codes
rows2 = (strfind(reshape(incident_symptom_code(:,2:3)',[1 length(incident_symptom_code)*2]),' (')+1)/2;
system_label = unique(Incident.Content(rows2,5));
tmp = char(system_label);
system_opt = system_label{strfind(tmp(:,1)',opt)}
filename = [system_opt '_Incidents_' datestr(Incident.start_date) '_to_' datestr(Incident.end_date) '.xlsx'];
xlswrite([folder '\' filename],Incident.Fieldname)
xlswrite([folder '\' filename],Incident.Content(rows,:),1,'A2')