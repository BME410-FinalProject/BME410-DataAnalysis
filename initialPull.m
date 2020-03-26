clc
clear

%%Import relevant data
amp = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\amplitudes.npy'); %Spike amplitude
time = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\spike_times.npy'); %Spike times
spiketemp = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\spike_templates.npy');
chanMap = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\channel_map.npy');%sorts channel with cluster
cluster = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\spike_clusters.npy');%spikes per cluster?

%%Plotting
figure(1)
plot(time(1:1000), amp(1:1000));