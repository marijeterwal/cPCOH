% computeConsensus

% Description:
% This function is called by computeCPCOH to compute the consensus PCOH
% after all permutations of the consensus-based partial coherence have been run.
% NOTE: the code removes all temporary files (see last section).

% INPUT: 
% The code requires two input arguments:
    % cfg: a config struct specifying the main parameters of the consensus analysis
    % (see ExampleScript.m for usage and checkConfig.m for defaults);
    % data: a FieldTrip data struct with preprocessed data.

% OUTPUT: 
% The resulting COH and cPCOH for each pair are stored in a data struct that is
% saved to a .mat file in the folder [cfg.path]/PartialCoherence/.

% This code implements the consensus-based partial coherence method as
% described in Ter Wal et al., NeuroImage, 2018.
% DOI: https://doi.org/10.1016/j.neuroimage.2018.06.011
% For more details and citations please refer to the paper. 

% Marije ter Wal, 2018
% m.terwal@donders.ru.nl || m.j.terwal@bham.ac.uk

function computeConsensus(cfg, data)

if isempty(cfg); error('Config missing'); end
if nargin < 2 || isempty(data); error('Data missing'); end

%% check config

cfg         = checkConfig(cfg, data, 'PCOH');
nevents     = length(cfg.events);
nperms      = cfg.nperms;

% channels & pairs
pairs = cfg.pairs;
npairs = size(pairs,1);

Zthres = abs(norminv(cfg.alpha/2,0,1));

%% permutations

for pr = 1:npairs
    
    fprintf('\nPair %i of %i... \n', pr, npairs)
    
    if ~cfg.overwrite && exist([cfg.path, 'PartialCoherence/', cfg.subject,...
            cfg.fileName, ...
            '_Channel', data.label{pairs(pr,1)}, '_Channel', data.label{pairs(pr,2)}, ...
            '_dataPCoh', '.mat'], 'file')
        continue
    end
    
    % preallocate
    dataPCoh = createCPCOHstruct(cfg,data,pairs(pr,:),1:nperms);
    
    for pm = 1:nperms
        
        %% save PCoh data per pair
        dum = loadpar([cfg.path, 'PartialCoherenceTempFiles/', cfg.subject,...
            cfg.fileName,...
            '_Channel', data.label{pairs(pr,1)}, '_Channel', data.label{pairs(pr,2)}, ...
            '_perm#' num2str(pm),'.mat']);
        permdum = dum.dataPCoh;
        
        %% transfer data to struct
        for ev = 1:nevents
            dataPCoh.Coh(ev,pm,:,:) = permdum.Coh(ev,:,:,:);
            dataPCoh.PCoh(ev,pm,:,:) = permdum.PCoh(ev,:,:,:);
            dataPCoh.CohZ(ev,pm,:,:) = permdum.CohZ(ev,:,:,:);
            dataPCoh.PCohZ(ev,pm,:,:) = permdum.PCohZ(ev,:,:,:);
        end
    end
    
    %% compute consensus
    for ev = 1:nevents
        datZ = dataPCoh.PCohZ(ev,:,:,:);
        datsign = zeros(size(datZ));
        datsign(abs(datZ) >= Zthres) = 1;
        dataPCoh.Consensus(ev,1,:,:) = nanmean(datsign,2);
        dataPCoh.PCohc(ev,1,:,:) = nanmean(dataPCoh.PCoh(ev,:,:,:),2) .* ...
            cfg.weightFunction(dataPCoh.Consensus(ev,:,:,:));
        dataPCoh.PCohZc(ev,1,:,:) = nanmean(dataPCoh.PCohZ(ev,:,:,:),2) .* ...
            cfg.weightFunction(dataPCoh.Consensus(ev,:,:,:));
    end
    
    %% save PCoh data per pair
    savepar([cfg.path, 'PartialCoherence/', cfg.subject,...
        cfg.fileName, ...
        '_Channel', data.label{pairs(pr,1)}, '_Channel', data.label{pairs(pr,2)}, ...
        '_dataPCoh.mat'], dataPCoh,'dataPCoh');
    
    %% test and delete temp files
    
    try
        % test
        loadpar([cfg.path, 'PartialCoherence/', cfg.subject,...
            cfg.fileName, ...
            '_Channel', data.label{pairs(pr,1)}, '_Channel', data.label{pairs(pr,2)}, ...
            '_dataPCoh.mat']);
    catch
        warning('Loading check for pair %i failed', pr)
        continue
    end
    
    % delete temp files
    if cfg.deleteTemp
        fprintf('Deleting temp files...\n')
        for pm = 1:nperms
            delete([cfg.path, 'PartialCoherenceTempFiles/', cfg.subject,...
                cfg.fileName, ...
                '_Channel', data.label{pairs(pr,1)}, '_Channel', data.label{pairs(pr,2)}, ...
                '_perm#' num2str(pm),'.mat']);
        end
    end

end % pairs

end

