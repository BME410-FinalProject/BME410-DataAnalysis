%% Calculate neural response (# spikes) to visual stimulus presentations
% Load sample data
% data downloaded from: http://data.cortexlab.net/singlePhase3/

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

figure, imagesc(resp)
title('Response Plot')

%%Z-Scored Data

z_resp = zscore(resp);

% Plot transposed matrix (makes it easier to see) and compare
figure
set(gcf,'Position',[100 500 1000 800])
imagesc(z_resp')
colorbar
title('Matlab z-score for Response Plot')
xlabel('Stimulus #')
ylabel('Neuron #')


%% We don't have working tuning curve data. Not sure if that is something
% we even need??


%% PCA on response data


[coeff,score,latent] = pca(resp);

figure
hold on
for k = 1:ncellsub
    scatter(score(1:ncellsub,1),score(1:ncellsub,2),'.');
end
xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')
hold off

% Scrap Code for testing purposes
% figure
% hold on
% % for k = 1:ncellsub
% scatter3(score(:,1),score(:,2),score(:,3),'.');
% % end
% xlabel('PC 1')
% ylabel('PC 2')
% zlabel('PC 3')
% hold off

figure
hold on
for k = 1:ncellsub
    scatter3(score(1:ncellsub,1),score(1:ncellsub,2),score(1:ncellsub,3),'.');
end

xlabel('PC 1')
ylabel('PC 2')
zlabel('PC 3')

%% Linear Regression on PCA Data
subN = round(length(clustids)/10); %%This is a guess
[b,bint,r,rint,stats] = regress(StimTimes', [ones(size(score,1),1),score(:,1:subN)]);

% OK, given our estimate of B, we can estimate stimuli variables 
% on single trials (y_i) based on the neural responses xi and b
% y_hat_i = bhat0 + bhat1*x1_i + bhat2*x2_i + ..... + bhatn*xn_i
statshat = [ones(size(score,1),1),score(:,1:subN)]*b; % y_hat = X*bhat;

figure
plot(StimTimes',statshat,'k.')
hold on
plot(StimTimes',StimTimes','r-')
legend('Estimates','unity')
xlabel('True Stim ID')
ylabel('Estimated Stim ID')
title(['R^2 =',num2str(stats(1)),' using ', num2str(subN), ' PCs']); % stats returns the R2 statistic, 
%the F-statistic and its p-value, and an estimate of the error variance.

%% Now, just use first 3 PCs
[b3,bint3,r3,rint3,stats3] = regress(StimTimes', [ones(size(score,1),1),score(:,1:3)]);
statshat3 = [ones(size(score,1),1),score(:,1:3)]*b3; % y_hat = Xb;
figure
plot(StimTimes',statshat3,'k.')
hold on
plot(StimTimes',StimTimes','r-')
legend('Estimates','unity')
xlabel('True Stim ID')
ylabel('Estimated Stim ID')
title(['R^2 =',num2str(stats3(1)),' using 3 PCs']); % stats returns the R2 statistic, 
%the F-statistic and its p-value, and an estimate of the error variance.

%% How does changing the number of PCs change the fit?
% choose a maximum number of PCs to be 1/3 of total trials
MaxN = round(length(StimTimes')/3);
% initialize a vector to store the R_sq for data as a function of # PCs
statschange = NaN(MaxN,1);
% go through and do linear regression to estimate B with a different number
% of PCs
MaxN = 50; %% I had to change this because it gets too good too fast
for k = 1:MaxN
    [b,bint,r,rint,stats] = regress(StimTimes', [ones(size(score,1),1),score(:,1:k)]);
    statschange(k) = stats(1);
end
figure
plot(1:MaxN,statschange(1:50))
xlabel('# PCs used in Regression')
ylabel('R^2 value of fit')
title('How does the # of PCs used in Regression affect prediction quality?')

%% To do "offline" neural decoding, we can train on a subset of data, 
% and "decode" the remainder
% initialize a vector of "training fractions"
fract_vect = [0.5:0.01:0.99];
% RSS = residual sum of squares = sum_i [(estimate - true)^2]
Rss_decode = NaN(length(fract_vect),1);
Ntrials = length(StimTimes);

% loop through each training fraction
for k = 1:length(fract_vect)
    % calculate # of training trials
Ntrain = round(fract_vect(k)*Ntrials);
% testing trials are the remainder
Ntest = Ntrials - Ntrain;

% Train our linear regression decoder
% training y = y data for training trials
% training x = neural responses on training trials (plus column of ones)
[btr,bint,r,rint,statstr] = regress(StimTimes(1:Ntrain)', [ones(Ntrain,1),score(1:Ntrain,1:MaxN)]);

% Test on left-out data
Xtest = [ones(Ntest,1),score(Ntrain+1:Ntrials,1:MaxN)];
ytest = StimTimes(Ntrain+1:Ntrials)';
ytesthat = Xtest*btr;
Rss_decode(k) = sum((ytest - ytesthat).^2); % compute residual sum squares
end

figure
plot(fract_vect,Rss_decode,'.-')
xlabel('Fraction of data used for training')
ylabel('Residual sum of squares')
title('Prediction errors as a function of training data')

%% Let's go back to thinkg about doing regression just on Neurons? 
% How should we pick them? Maybe choose those with the highest FR?
% To begin with, we ignore tuning 

avgFR = mean(resp);
[sortedFR,sortidx] = sort(avgFR,'descend');
% switch response matrix acorrding to sort
sortedNeur = resp(:,sortidx); 
% neural noise is not gaussian, BUT
% linear regression assumes gaussian noise; taking the square root of the 
% responses makes them look more gaussian
sortedNeur_sqrt = sqrt(sortedNeur); 

% Not sure about the tuning parts
% Now choose TUNED neurons with highest FR
% [sortedtuned,sortidxtune] = sort(mean(resptuned),'descend');
% sortedNeurtuned = resptuned(:,sortidxtune);
% sortedNeurtunedsqrt = sqrt(sortedNeurtuned);

% what if look at zscored data w/o pca?
[sortedz,sortidxz] = sort(mean(z_resp));
sortedztuned = z_resp(:,sortidxz);

% Lets run linear regression on Neuron data
% How does changing the number of PCs change the fit?
%MaxN = round(length(stim.istim)/3); % same maxN as for PC data
MaxN = 25;
statschangeNEU = NaN(MaxN,1);
statschangeNEUsqrt = NaN(MaxN,1);
statsNeurtune = NaN(MaxN,1);
statsNeurtunesqrt = NaN(MaxN,1);
statsz = NaN(MaxN,1);
for k = 1:MaxN
    [~,~,~,~,stats] = regress(StimTimes', [ones(size(sortedNeur,1),1),sortedNeur(:,1:k)]);
    statschangeNEU(k) = stats(1);
    [~,~,~,~,statssqrt] = regress(StimTimes', [ones(size(sortedNeur_sqrt,1),1),sortedNeur_sqrt(:,1:k)]);
    statschangeNEUsqrt(k) = statssqrt(1);
    %[~,~,~,~,statstuned] = regress(StimTimes', [ones(size(sortedNeurtuned,1),1),sortedNeurtuned(:,1:k)]);
    %statsNeurtune(k) = statstuned(1);
    %[~,~,~,~,statstunedsqrt] = regress(StimTimes', [ones(size(sortedNeurtunedsqrt,1),1),sortedNeurtunedsqrt(:,1:k)]);
    %statsNeurtunesqrt(k) = statstunedsqrt(1);
    [~,~,~,~,statszz] = regress(StimTimes', [ones(size(sortedztuned,1),1),sortedztuned(:,1:k)]);
    statsz(k) = statszz(1);
end

figure
plot(1:MaxN,statschange(1:MaxN))
hold on
plot(1:MaxN,statschangeNEU)
plot(1:MaxN,statschangeNEUsqrt)
plot(1:MaxN,statsNeurtune)
plot(1:MaxN,statsNeurtunesqrt)
plot(1:MaxN,statsz)
legend('PCs','Neurons','sqrt(Neu)','tuned Neurons','tuned Neurons sqrt','z-scored data')
xlabel('# PCs/Neurons used in Regression')
ylabel('R^2 value of fit')
title('PCA vs. single neuron response amplitudes?')