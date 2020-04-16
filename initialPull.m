clc
clear

%%Import relevant data
%{
Note that in order to import data, you will need to change the directory
to refer to where you have the data stored on your personal computer. I
don't know if there is a better way to do this, but here's a template at
the very least. It's easiest if you navigate to the file using your file
explorer, and then copy the directive from the top bar. You have to add in
the file name manually with this method
%}
amp = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\amplitudes.npy'); %Spike amplitude
raw_time = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\spike_times.npy'); %Spike times (by sample)
chanMap = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\channel_map.npy');%sorts channel with cluster
cluster = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\spike_clusters.npy');%spikes per cluster?
chanPos = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\channel_positions.npy');

%% Stimulus
load('stimInfo.mat');
stimulusTime = stimTimes{1,1};
stimType = stimPositions{1,1};

%% Manipulate imported data
time_ms = double(raw_time)/30000.0; % Divide by sampling rate to get in milliseconds
numClust = max(cluster); %number of clusters ("neurons") determined
chan_split = zeros(numClust, size(stimulusTime, 2)); 
time_split = zeros(numClust, size(stimulusTime, 2)); 

%% Separate spikes per cluster for plotting
for i = 1:size(stimulusTime, 2)
    chan_split(cluster(i)+1, i) = amp(i);
    time_split(cluster(i)+1, i) = time_ms(i);
end

%% Plotting spike data
figure;
plot(time_split(88, :), chan_split(88, :), "o"); %plotting spikes for cluster 87 (88 minus 1)
title('Spike amplitudes over time');
xlabel('Time (s)');
ylabel('Amplitude (mV?)');

%% Z score data
z_scored_data = zscore(time_split);

figure;
subplot(1, 2, 2);
imagesc(z_scored_data);
colorbar;
title('Z-Scored spike time data for 297 spike clusters')
xlabel('Time')
ylabel('Cluster number')

subplot(1, 2, 1);
imagesc(time_split(:,:));
colorbar;
title('Raw spike time data for 297 spike clusters');
xlabel('Time (ms)');
ylabel('Cluster number');

%% "Tuning Curve"
coor = unique(stimType, 'rows'); %Number of unique coordinates tested
Ncoor = length(coor); %How many different stimuli?
row = max(coor(:, 1));
col = max(coor(:, 2));
tuning = NaN(numClust, row, col);

nums = round(1 + (297-1) .* rand(5,1)); %generate random neurons in our range to investigate

for i = 1:5 %change to Nneurons if you want all neurons
    temp = zeros(row, col); % reset temporary array
    for j = 1:size(time_split, 2) %iterate through length of time
        if (time_split(nums(i), j) > 0) %check whether or not there is a spike at this time
            k = stimType(j, 1); %check which stimulus position triggered the spike
            m = stimType(j, 2);
            temp(k, m) = temp(k, m) + 1;
        end
    end 
    for k = 1:row
        for m = 1:col
            tuning(nums(i), k, m) = temp(k, m); %Fill tuning array accordingly
        end
    end
end

figure; 
hold on
tempPlot = zeros(row, col);
count = 1;
for i = 1:5
    for k = 1:row
        for m = 1:col
            tempPlot(k, m) = tuning(nums(i), k, m);
        end
    end
    subplot(2, 3, count);
    imagesc(tempPlot);
    colorbar;
    neuron = nums(i);
    title(['Plot of stimulus response for neuron ' int2str(neuron)]);
    count = count + 1;
end
hold off