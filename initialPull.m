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
time = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\spike_times.npy'); %Spike times (NOT SORTED LINEARLY)
spiketemp = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\spike_templates.npy');
chanMap = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\channel_map.npy');%sorts channel with cluster
cluster = readNPY('C:\Users\karab\Documents\GitHub\BME410-DataAnalysis\spike_clusters.npy');%spikes per cluster?

%%Plotting
figure(1)
plot(time(1:1000), amp(1:1000)); %only plotting the first 1000 datapoints to see shape