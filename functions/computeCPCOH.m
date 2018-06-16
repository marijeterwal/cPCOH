
% computeCPCOH

% Description:
% This function takes a FieldTrip data struct and computes the coherence
% and consensus-based partial coherence for pairs and channels indicated in
% the cfg.

% INPUT: 
% The code requires two input arguments:
    % cfg: a config struct specifying the main parameters, such as the number 
    % of groups and permutations (see ExampleScript.m for usage and checkConfig.m 
    % for defaults);
    % data: a FieldTrip data struct with preprocessed data.

% The code allows for 3 different types of reference distibution for the 
% Z-scoring opperation: (text in brackets indicates how to choose each reference)
% - relative to baseline (specify cfg.baseline)
% - a trial shuffle of the samen data with another dataset (specify cfg.nreps)
% - a different reference dataset, for example a rest condition or
% independently recorded baseline (specify neither cfg.baseline nor cfg.nreps)

% OUTPUT: 
% The resulting COH and cPCOH for each pair are stored in a data struct that is
% saved to a .mat file in the folder [cfg.path]/PartialCoherence/. Intermediate 
% results are stored in the folder [cfg.path]/PartialCoherenceTempFiles/, 
% to allow for restarts without losing work. 
% When only 1 pair is computed, the data struct will be returned. 


% This code implements the consensus-based partial coherence method as
% described in Ter Wal et al., NeuroImage, 2018.
% DOI: https://doi.org/10.1016/j.neuroimage.2018.06.011
% For more details and citations please refer to the paper. 

% Marije ter Wal, 2018
% m.terwal@donders.ru.nl || m.j.terwal@bham.ac.uk


function dataOut = computeCPCOH(cfg, data, dataRef)

%% check paths
if ~exist([cfg.path, 'PartialCoherence/'], 'dir')
    dirs = input(['PCOH directories required for saving the results were not found,'...
        'would you like to create them? (y/n)\n'], 's');
    if strcmp(dirs,'y')
        mkdir([cfg.path, 'PartialCoherence/'])
        mkdir([cfg.path, 'PartialCoherenceTempFiles/'])
    else
        fprintf('Missing directories, cPCOH computation aborted\n')
        return
    end
end

%% check config and find Zs-scoring reference
if nargin == 2; dataRef = []; end
if nargin < 2 || isempty(data); error('Data missing'); end

% check config
cfg         = checkConfig(cfg, data, 'PCOH');

% determine type of reference
if isfield(cfg, 'baseline') % baseline in the same dataset
    refType = 'Baseline';
    fprintf('Baseline reference selected\n')
    cfg.nreps = length(cfg.baseline(1):cfg.dt:cfg.baseline(2));
elseif isfield(cfg, 'nreps')
    refType = 'TrialShuffle'; % trial shuffling
    fprintf('Trial shuffle reference selected\n')
else
    refType = 'Reference';  % reference dataset
    fprintf('Reference dataset selected\n')
    cfg.nreps = size(dataRef.time{1});
end

if (strcmp(refType, 'TrialShuffle') || strcmp(refType, 'Reference')) && isempty(dataRef)
    error('No reference dataset specified')
end

%% set some general values
nfreqs      = length(cfg.foi);
ntimes      = length(cfg.toi);
nevents     = length(cfg.events);
ntrials     = length(data.trial);
nreps        = cfg.nreps;

%% channels and pairs
channels = unique([cfg.channels(:);cfg.pairs(:)]);
pairs = cfg.pairs;

[~, pairsID(:,1)] = ismember(pairs(:,1), channels);
[~, pairsID(:,2)] = ismember(pairs(:,2), channels);

npairs = size(pairs,1);
nchannels = length(channels);

fprintf('Analyzing %i pairs\n', npairs);

if nchannels < max(cfg.ngroups); error('Not enough channels for the number of groups'); end

%% define conditioning sets
CondGroups = getPermutations(cfg);
cfg.nperms = length(CondGroups);
nperms = cfg.nperms;

%% compute spectra

fprintf('Test dataset: Computing spectra...')
cfs1c = getSpectra(cfg, data);
fprintf(' Done!\n')

fprintf('Reference dataset: Computing spectra...')
switch refType
    case 'TrialShuffle'
        cfs2c = getSpectra(cfg, dataRef);
    case 'Reference'
        cfgt = cfg;
        cfgt.events = 1;
        cfs2c = getSpectra(cfgt, dataRef);
    case 'Baseline'
        cfgt = cfg;
        cfgt.events = 1;
        cfgt.toi = cfg.baseline(1):cfg.dt:cfg.baseline(2);
        cfs2c = getSpectra(cfgt, data);
end
fprintf(' Done!\n\n')

%% permutations

condGr = {CondGroups.Groups};
saveQual = ones(nperms,1);

% This code is compatible with the Parallel Computing Toolbox: 
% to run the code in parallel change the outer for-loop to parfor
for pm = 1:nperms
    fprintf('\nPermutation %i of %i... ', pm, nperms)
    
    for pr = 1:npairs
        
        % does pcoh already exist?
        if ~cfg.overwrite
            if exist([cfg.path, 'PartialCoherenceTempFiles/', cfg.subject,...
                    cfg.fileName, ...
                    '_Channel', data.label{pairs(pr,1)}, '_Channel', data.label{pairs(pr,2)}, ...
                    '_perm#' num2str(pm),'.mat'], 'file')
                continue
            end
            
            if exist([cfg.path, 'PartialCoherence/', cfg.subject,...
                    cfg.fileName, ...
                    '_Channel', data.label{pairs(pr,1)}, '_Channel',data.label{pairs(pr,2)},...
                    '_dataPCoh', '.mat'], 'file')
                continue
            end
        end
        
        % preallocate
        dataPCoh = createCPCOHstruct(cfg,data,pairs(pr,:),pm);
        
        for nf = 1:nfreqs
            % slice variables
            cfs1 = cfs1c(:,:,:,nf,:);
            cfs2 = cfs2c(:,:,:,nf,:);
            
            %---------- calculate Partial Coherence
            % determine partial coherence
            for ev = 1:nevents
                
                % cfs: trials x events x channels x freqs x times
                [Coh,PCoh,Qtk] = spectra2CPCOH(squeeze(cfs1(ev,:,:,:,:)),pairsID(pr,:), condGr{pm});
                
                % if conditioning doesn't work, move on to next conditioning permutation
                if Qtk == 0; saveQual(pm) = 0; warning('Bad permutation'); break; end
                
                % store data
                dataPCoh.PCoh(ev,1,nf,:) = PCoh; % save trial data
                dataPCoh.Coh(ev,1,nf,:) = Coh;
            end
            if saveQual(pm) == 0; break; end % break trial loop if conditioning failed
            
            %---------- calculate Reference
            if nreps < 5 % don't do shuffling is the number of repetitions is too low for statistics
                warning('Z-scoring was skipped due to low number of repetitions')
                continue
            end
            
            switch refType
                case 'TrialShuffle'
                    eventRef = 1:nevents;
                    [Coh0Ev,PCoh0Ev] = deal(zeros(nevents,nreps, ntimes));
                    
                    % randomly choose trials
                    for rp = 1:nreps
                        % make randomized dataset
                        cfsrp = cfs2;
                        for nc = 1:nchannels
                            if randn < 0.5
                                cfsrp(:,:,nc,:,:) = cfs1(:,randperm(ntrials),nc,:,:);
                            else
                                cfsrp(:,:,nc,:,:) = cfs2(:,randi(length(dataRef.trial),[ntrials,1]),nc,:,:);
                            end
                        end
                        cfsrp(:,:,pairsID(pr,1),:,:) = cfs1(:,randperm(ntrials),pairsID(pr,1),:,:);
                        cfsrp(:,:,pairsID(pr,2),:,:) = cfs2(:,randi(length(dataRef.trial),[ntrials,1]),pairsID(pr,2),:,:);
                        
                        % partial coherence on shuffled data
                        for ev = 1:nevents
                            [Coh0,PCoh0,Qtk] = spectra2CPCOH(squeeze(cfsrp(ev,:,:,:,:)),pairsID(pr,:), condGr{pm});
                            if Qtk == 0; saveQual(pm) = 0; break; end
                            
                            PCoh0Ev(ev,rp,:) = PCoh0;
                            Coh0Ev(ev,rp,:) = Coh0;
                        end
                        if saveQual(pm)== 0; warning('Bad permutation'); break; end
                    end
   
                case {'Reference', 'Baseline'}
                    eventRef = ones(nevents,1);
                    [Coh0,PCoh0,~] = spectra2CPCOH(squeeze(cfs2(1,:,:,:,:)),pairsID(pr,:), condGr{pm});
                    if saveQual(pm)== 0; warning('Bad permutation'); continue; end
                    
                    PCoh0Ev = reshape(PCoh0, [1,nreps,1]);
                    Coh0Ev = reshape(Coh0, [1,nreps,1]);
            end
            
            %---------- calculate z-score - Rank sum test
            for ev = 1:nevents
                dataPCoh.CohZ(ev,1,nf,:) = computeZscore(dataPCoh.Coh(ev,1,nf,:), Coh0Ev(eventRef(ev),:,:));
                dataPCoh.PCohZ(ev,1,nf,:) = computeZscore(dataPCoh.PCoh(ev,1,nf,:), PCoh0Ev(eventRef(ev),:,:));
            end
            
        end % freqs
        
        
        %% save PCoh data per pair
        savepar([cfg.path, 'PartialCoherenceTempFiles/', cfg.subject,...
            cfg.fileName, ...
            '_Channel', data.label{pairs(pr,1)}, '_Channel', data.label{pairs(pr,2)}, ...
            '_perm#' num2str(pm),'.mat'], dataPCoh, 'dataPCoh');
        
    end % pairs
    
end % permutations

% delete(pp)

%% find consensus and merge files into one output file
computeConsensus(cfg, data)

%% output data

% store quality check of permutations
dum = num2cell(saveQual);
[CondGroups.Quality] = dum{:};
save([cfg.path, 'PartialCoherence/', cfg.subject, ...
        cfg.fileName, '_CondGroups','.mat'], 'CondGroups');

% return or not
if npairs == 1
    dum = loadpar([cfg.path, 'PartialCoherence/', cfg.subject,...
            cfg.fileName, ...
            '_Channel', data.label{pairs(1,1)}, '_Channel', data.label{pairs(1,2)}, ...
            '_dataPCoh.mat']);
        
    dataOut = dum.dataPCoh;
else
    warning('Too many pairs to output; load individual pairs instead.')
    dataOut = NaN;
end

end



function [Z] = computeZscore(dat, ref)

% dat = 1 x 1 x 1 x ntimes
% ref = 1 x nrep x ntimes or ref = 1 x nrep x 1

if size(ref,4)==1 && size(ref,4)>1
    error('Ref data is of the wrong size')
elseif size(ref,3)>1 && size(dat,4) ~= size(ref,3)
    error('Data and Ref data are of different length')
end

ntimes = size(dat,4);
Z = squeeze((dat - repmat(mean(ref,2), [1,1,1,ntimes])) ./ repmat(std(ref,1,2),[1,1,1,ntimes])); 

end

