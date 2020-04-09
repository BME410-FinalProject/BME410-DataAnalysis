%{
BME 410 Data Analysis Project
Code Submission 4.2.20
Noah Smith, Charlie Martin, Lindsay Karaba

This program extracts the raw voltage data recorded from neuropixel
probes and forms two plots. The first shows voltage (represented by color)
across 72 channels over the timespan. The second extracts the voltage trace
over time for selected channels of the probe.
%}
%% Ensure numpy path is established
%note that addpath path must be changed for specific computer
addpath('D:\Documents\Word\College\Junior Year\Semester 2\BME 410\npy-matlab-master')
savepath

%% Read in raw data file
%note that fn path must be changed for specific computer
fn = 'D:\Documents\Word\College\Junior Year\Semester 2\BME 410\rawDataSample.bin';
fid = fopen(fn, 'r');
rawData = fread(fid, [385 Inf], '*int16');
fclose(fid);
chanMap = readNPY('channel_map.npy');
rawData = rawData(chanMap+1,:);

%% Plot color map of data
figure;
imagesc(rawData(:,:));
colorbar;
title('Raw voltage data for 72 channels over time');
xlabel('Time (ms)');
ylabel('Channel number');

%% Plot set of traces
figure;
hold on
for i = 1:5 %arbitrary number of traces, chose 5 for simplicity (max of 72)
    subplot(5, 1, i);
    plot(1705000:1800000, rawData(i, 1705000:1800000)); %plot small section of time
    if i < 5
        set(gca, 'xticklabel', []) %remove xticks for all but bottom plot
    end
    if (i == 1)
        title('Raw voltage traces of channels 1-5 (descending) over time');
    end
end
xlabel('Time (ms)')
ylabel('Voltage (mV)')
hold off

%% Now need to add in the Z-scored data
raw_double = double(rawData);
z_scored_data = zscore(raw_double);


figure;
subplot(1, 2, 2);
imagesc(z_scored_data);
colorbar;
title('Z-Scored Voltage Data for 72 channels')
xlabel('Time')
ylabel('Channel Number')

subplot(1, 2, 1);
imagesc(rawData(:,:));
colorbar;
title('Raw voltage data for 72 channels over time');
xlabel('Time (ms)');
ylabel('Channel number');

%% Import stimulus information
load('stimInfo.mat');
stimulusTime = stimTimes{1,1};
stimType = stimPositions{1,1};
response = rawData(:, size(stimulusTime, 2));
z_scored_short = rawData(:, size(stimulusTime,2));

%% Tuning
coor = unique(stimType, 'rows'); %Number of unique coordinates tested
Ncoor = length(coor); %How many different stimuli?
Nneurons = size(response,1); %How many neurons?
row = size(coor, 1);
col = size(coor, 2);
temp = NaN(size(response,2));
m = 1;
hold on
figure;
for i = 1:row
    for j = 1:col
        for k = 1:size(stimulusTime, 2)
           subplot(row, col, m);
           if (stimType(k, 1) == (i))&&(stimType(k, 2) == (j)) 
               temp(k) = response(k);
           end
           tempPlot = mean(temp(:));
           imagesc(stimulusTime, tempPlot);
           colorbar;
        end
    end
end
           
               
% tuningcurves = NaN(Ndir,Nneurons); 
% ztuning = NaN(Ndir,Nneurons); 
% for k = 1:Ndir % go throuh each direction    
%     tuningcurves(k,:) = mean(response(stimType == dir(k, :),:)); % calculate mean response of each neuron to single stim direction    
%     ztuning(k,:) = mean(z_scored_short(stimType == dir(k, :),:)); % "" for z-scored data 
% end
% figure(3) 
% set(gcf,'Position',[100 500 1000 500]) 
% subplot(2,1,1) 
% imagesc(tuningcurves) 
% xlabel('Neuron #') 
% ylabel('Mean resp.') 
% title('Tuning Curves') 
% subplot(2,1,2)
% imagesc(ztuning) 
% xlabel('Neuron #') 
% ylabel('Mean resp.') 
% title('Tuning Curves of Z-scored data')
% % claim neurons are tuned if at any orientation have z-score > 0.5
% [maxztc,PD] = max(ztuning); 
% tuned = find(maxztc >= 0.5);

