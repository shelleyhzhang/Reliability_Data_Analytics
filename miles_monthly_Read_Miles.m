function Read_Miles_monthly
% Author: Shelley Zhang
% Email:  szhang@ltk.com
% Date:   Sep 2017
%%
clear all
close all
%To manually select folder
%folder = uigetdir('P:\C4314 MBTA-X29PS01 Evaluation of Light Rail Vehicles\Type 8 Truck Overhaul\MCRS Analysis\Type 8_Mileage by Unit\');
folder = 'P:\C4314 MBTA-X29PS01 Evaluation of Light Rail Vehicles\Type 8 Reliability Analysis\MCRS\GLT8_Mileages Monthly\';
folder = 'C:\Users\SZhang\Documents\MATLAB\T8\GLT8_Mileage_Unit_Monthly\';
filelist = ls([folder '\20*.xls*']);
date_first = '2015-07-01';
date_last = '2018-12-31'; 
Intervals = size(filelist,1);
%%%%%%%%%% Prior Knowledge of the Fleet T8%%%%%%%%%%
Mileage.Unit = [3800:3894];
Mileage.Vals = zeros(length(Mileage.Unit),Intervals);
Mileage.start_date = datenum(date_first);
Mileage.end_date = datenum(date_last);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n=1:Intervals
    [num, text] = xlsread([folder '\' filelist(n,:)]);
    %in individual nth month
    Unitn = str2num(cell2mat(text(2:end-1,1)));
    index = find(ismember(Mileage.Unit, Unitn)); %note Car.Unit is larger set than Unitn
    Mileage.Vals(index,n) = num(2:end,1);
    fprintf('Read %2g:%s\n',n, filelist(n,:))
end
save([folder '\..\DATA_mat\Mileages_' date_first '_to_' date_last '.mat'],'Mileage')
%% To Print Figure
load([folder '\..\DATA_mat\Mileages_' date_first '_to_' date_last '.mat'],'Mileage')
% tmp = findstr(filelist(1,:),'_Q');
% start_time = filelist(1,tmp-4:tmp+2);
first_unit = Mileage.Unit(1);
figure;
imagesc(Mileage.Vals),colorbar
xticks(1:Intervals);
% tmp = datestr((ones(Intervals,1)*datenum('04/01/2010')+[0:92:92*(Intervals-1)]'));
% tmp = datetime(2010,4:4:Intervals*4+4,1)
category_in_quarter = discretize(datetime(2010,4:3:(Intervals*3+4),1),'quarter','categorical');
set(gca,'XTickLabel',char(category_in_quarter),'XTickLabelRotation',90)
%addnote(start_time,first_unit,Intervals);

%to create an array of datetime from Q2 2010, should change if start time is different

%to mask any travel miles less than 50 quarterly in red
tmp = find(Mileage.Vals<50);
mask = zeros(size(Mileage.Vals)); mask(tmp)=1;
figure,imagesc(cat(3,Mileage.Vals/max(Mileage.Vals(:))+mask,Mileage.Vals/max(Mileage.Vals(:)),Mileage.Vals/max(Mileage.Vals(:))))
addnote(date_first,first_unit,Intervals);
xticklabels(char(category_in_quarter))

% to generate histogram of qualified quarterly miles for all cars
figure,h=histogram(Mileage.Vals(Mileage.Vals>50),100);xlabel('MMBF [Mile]'),ylabel('Occurences [Quarter]')
%h = histogram(Car.Mile(Car.Mile>50),100);
count = h.Values;
center = h.BinEdges(1:end-1)+diff(h.BinEdges)/2;
hold on
tmp = get(gca,'xlim'); y = linspace(tmp(1),tmp(2),1000); 
[tmp tmp2] = sort(count); mu = mean(center(tmp2(end-2:end)));
[tmp] = find(count>0.5*max(count)); sigma = (max(center(tmp))-min(center(tmp)))/2.355; %FWHM is 2.355 times of sigma
f = max(count)*exp(-(y-mu).^2./(2*sigma^2));
plot(y,f,'k--','LineWidth',1.5)
%hold on, plot(center, count)
function addnote(start_time,first_unit,Intervals)
title('MMBF by Car (-Y) by Quarter (-X) [Mile]')
%xlabel(['Intervals starting from ' start_time ' [Quarter]'],'Interpreter','none')
xlabel(['Intervals [Quarter]'])
set(gca,'XTickLabelRotation',45)
ylabel('Car [Unit]'),colormap(gray)
tmp = get(gca,'YTick')+first_unit-1; 
set(gca,'YTickLabel',tmp,'fontsize',18,'fontweight','bold')
tmp = [1:Intervals];set(gca,'XTick',[1:Intervals]), set(gca,'XTickLabel',tmp,'fontsize',18,'fontweight','bold')
colormap(specialcmap('rgc'))
