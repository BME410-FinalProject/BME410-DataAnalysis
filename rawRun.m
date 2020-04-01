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
