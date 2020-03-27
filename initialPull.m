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
raw_time = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\spike_times.npy'); %Spike times (NOT SORTED LINEARLY)
spiketemp = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\spike_templates.npy');
chanMap = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\channel_map.npy');%sorts channel with cluster
cluster = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\spike_clusters.npy');%spikes per cluster?
chanPos = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\channel_positions.npy');

%%Manipulate imported data
time_s = raw_time./30000; % Divide by sampling rate to get in milliseconds
numClust = max(cluster); %number of clusters ("neurons") determined
chan_split = zeros(numClust, 198); 
time_split = zeros(numClust, 198); 

%% Separate spikes per cluster for plotting
for i = 1:58544
    chan_split(cluster(i)+1, i) = amp(i);
    time_split(cluster(i)+1, i) = time_s(i);
end

%%Plotting
figure(1)
plot(time_split(88, :), chan_split(88, :), "o"); %plotting spikes for cluster 87 (88-1)
title('Spike amplitudes over time');
xlabel('Time (s)');
ylabel('Amplitude (mV?)');
