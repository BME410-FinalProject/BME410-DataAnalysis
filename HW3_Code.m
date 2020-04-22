%% Calculate neural response (# spikes) to visual stimulus presentations
% Load sample data

% First, load stim info (stimulus positions and time of appearance)
% load('stimInfo.mat')
load('stimInfo.mat')
samp_rate = 30000;
StimPositions = stimPositions{1};
StimTimes = stimTimes{1}/samp_rate; % in (s), onset times
Nstim = length(StimTimes);

stim_duration = 1/6; % s

% Next, get spike times
samp_rate = 30000;
sptimes = double(readNPY('spike_times.npy'))/samp_rate;
%/samp_rate; 
% ^^ these times look like they are in terms of sample #, not actual time
% so convert to time by dividing by sample rate

% then, get ID of spikes (labels mapping spikes to single cells)
spclusters = readNPY('spike_clusters.npy');
clustids = unique(spclusters);
ncell = length(clustids);


% lets grab spike times for the first cell
sp = cell(ncell,1);
% sp = sptimes(spclusters == clustids(1));

% calculate # spikes for each stimulus (from on time to off time)
% Nsp = NaN(Nstim,1);
% for k = 1:Nstim
%     Nsp(k) = sum(sp >= StimTimes(k)  & sp < (StimTimes(k)+stim_duration));
% end
emptycells = NaN(ncell,1);
for k = 1:ncell
    sp{k} = sptimes(spclusters == k);
    emptycells(k) = isempty(sp{k});
end

% get rid of empty cells
spsub = sp(emptycells==0);
ncellsub = length(spsub);
% calculate # spike times

resp = NaN(Nstim,ncellsub);
for k = 1:ncellsub
    k;
    for g = 1:Nstim
        resp(g,k) = sum(spsub{k} > StimTimes(g) & spsub{k} < StimTimes(g)+stim_duration);
    end
end

figure(), imagesc(resp)
title('Response Plot')
xlabel('Neuron #')
ylabel('Time of Response')
colorbar

[coeff,score,latent] = pca(resp);

figure(2)
hold on
for k = 1:ncellsub
    scatter(score(1:ncellsub,1),score(1:ncellsub,2),'o');
end
title(sprintf('Percent of Variance Captured by PC1 and PC2 = %0.2f',sum(latent(1:2))/sum(latent)*100))
xlabel('PC 1')
ylabel('PC 2')
hold off

figure(3)
hold on
for k = 1:ncellsub
    scatter3(score(1:ncellsub,1),score(1:ncellsub,2),score(1:ncellsub,3),'o');
end
title(sprintf('Percent of Variance Captured by PC1, PC2, and PC 3 = %0.2f',sum(latent(1:3))/sum(latent)*100))
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')
hold off