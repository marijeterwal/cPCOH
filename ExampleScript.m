%%% Example script calling cPCOH

% Description:
% This script demonstrates the use of the computeCPCOH function and related 
% functions by computing the cPCOH for the model dataset ExampleData.mat.

% This code implements the consensus-based partial coherence method as
% described in Ter Wal et al., NeuroImage, 2018.
% DOI: https://doi.org/10.1016/j.neuroimage.2018.06.011
% For more details and citations please refer to the paper. 

% Marije ter Wal, 2018
% m.terwal@donders.ru.nl || m.j.terwal@bham.ac.uk

%% Clear workspace

% clear all
close all
clc

%% Set data path and load example data

path_in         =  'cPCOH\';
subject         = 'ExampleData6';
path_data       = [path_in, subject, '.mat'];

% Load data_Test
load(path_data)

%% Set config and compute Partial Coherence

% config
cfg             = [];
cfg.subject     = subject;
cfg.path        = path_in;  % specify if path for saving should be different from path_data
cfg.fileName    = '_Test';  % optional name for saving
cfg.overwrite   = true;     % if true, already existing analyses will be overwritten
cfg.deleteTemp  = true;     % if true, intermediate results will be deleted at completion

% specify reference:

% --- for using baseline as reference
cfg.baseline    = [-0.4, 0];  % s 

% --- for using a trial shuffling 
% (pass second dataset to computeCPCOH function)
% cfg.nreps       = 100;    % number of repetitions 

% --- for using a independent reference dataset
% (pass second dataset to computeCPCOH function)
% do not specify baseline and nreps

% specify cPCOH settings:
cfg.ngroups         = [2,6];     % number of groups: [min, max]
cfg.nperms          = 50;        % number of permutations
cfg.consensusThres  = 0.9;       % consensus threshold 

% weight function for computing masked PCOH 
% cfg.weightFunction  = @(x) max(0,sign(x-cfg.consensusThres)); % step function
cfg.weightFunction  = @(x) min(1,(1/cfg.consensusThres)*x); % linear decrease below threshold

cfg.wavelet         = 'cmor3-1'; % wavelet

% specify data to analyze:
cfg.pairs       = [[1,2]];       % specify pairs; emtpy [] for all pairs;
cfg.foi         = [30:2:70];  % frequencies of interest in Hz
cfg.toi         = -0.150:0.005:0.250;   % times of interest in s

% alternatively use:
% cfg.toilim      = [0, 0.3];%
% cfg.dt          = 0.005;  % s

% run cPCOH!
dataPCoh = computeCPCOH(cfg, data_Test);

%% Plot the results

figure('Position', [50,50,1000,500]); 
subplot(131)
imagesc(cfg.toi, cfg.foi, squeeze(dataPCoh.CohZ(1,1,:,:))); 
axis xy
colorbar; caxis([-3,3]);
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Coherence')

subplot(132)
imagesc(cfg.toi, cfg.foi, squeeze(abs(dataPCoh.Consensus(1,:,:,:)))); 
axis xy
colorbar; caxis([0,1]);
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Consensus')

subplot(133); hold on;
imagesc(cfg.toi, cfg.foi, squeeze(dataPCoh.PCohZc(1,:,:,:))); 
contour(cfg.toi, cfg.foi, abs(squeeze(dataPCoh.CohZ(1,1,:,:))), 1.9, 'LineColor', [1,1,1], 'lineWidth', 2);
axis xy
colorbar; caxis([-3,3]);
xlabel('Time (s)')
ylabel('Frequency (Hz)')
title('Masked PCOH Z-score')

