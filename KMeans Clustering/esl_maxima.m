%%% Extreme sea level 
%%% To 3-day maxima and save coordinates
%%% Outliers: station 7568, 8226

clear; close all; warning off; clc 

% add working folder (where data locates)

gen_folder= [pwd '\'];

data_folder = 'C:\Users\hli490\Desktop\Spatial_dependence\Input\data_gtsm_era5\cf_esl_daily_maxima\';

addpath(data_folder)

fsave = [gen_folder 'ESL_time_frames_and_zones\Updated\'];

zones = {'IO', 'NA', 'NEA', 'NWA', 'NWP', 'Oceania', 'SA', 'SEA', 'SWA', 'SWP'}; % Global zones


%% loading data 

zone = 'NA';
station_ids = xlsread(strcat('Zones\Updated\',zone,'_stations.xlsx'),'A:A');
station_num = length(station_ids);  
lonlat = zeros(station_num,2); % Coordinates

for i = 1 : station_num 
    station_name = sprintf('gtsm_station%05d.nc', station_ids(i));
    lonlat(i,1) = ncread(station_name, 'lon');
    lonlat(i,2) = ncread(station_name, 'lat');
    waterlevel(:,i) = ncread(station_name, 'waterlevel');   
end

time = ncread(station_name, 'index');

%% save

NA.station_ids = station_ids;
NA.time = time;
NA.waterlevel = waterlevel;
NA.lonlat = lonlat;

save([fsave 'NA.mat'],'NA','-mat','-v7.3','-nocompression');


