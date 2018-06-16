% computeConsensus

% Description:
% This function is called by computeCPCOH to generate a list of unique 
% conditioning sets (i.e. groups of channels to average over) for each of 
% the permutations of the consensus-based partial coherence.
% Firstly, a number is groups is randomly chosen from the values specified in
% cfg.ngroups. Secondly, the number of channels per group is generated
% randomly, with the restriction that groups cannot be empty.

% INPUT: 
% The code requires two input arguments:
    % cfg: a config struct specifying the main parameters, such as the number 
    % of groups and permuations (see ExampleScript.m for usage and checkConfig.m 
    % for defaults);
    % data: a FieldTrip data struct with preprocessed data.

% OUTPUT: 
% The functions outputs the struct CondGroups, which is also stored as a .mat
% file in [cfg.path]/PartialCoherence/ for future reference.

% This code implements the consensus-based partial coherence method as
% described in Ter Wal et al., NeuroImage, 2018.
% DOI: https://doi.org/10.1016/j.neuroimage.2018.06.011
% For more details and citations please refer to the paper. 

% Marije ter Wal, 2018
% m.terwal@donders.ru.nl || m.j.terwal@bham.ac.uk


function CondGroups = getPermutations(cfg)

%TODO: extend permutation set
% if CondGroups already exists, but nperm is larger than the number of permutations in CondGroups, then extend the existing set of groups to nperm and save under the old name 

try % load conditioning sets
    load([cfg.path, 'PartialCoherence/', cfg.subject, ...
        cfg.fileName, '_CondGroups','.mat']);
    fprintf('Found pre-existing conditioning set')
catch
    nchannels = length(cfg.channels);
    
    if max(cfg.ngroups) == length(cfg.channels)
        warning('Conditioning on all channels individually, npermutations set to 1.')
        CondGroups = struct('Permutation', num2cell(1), 'Ngroups', '', 'Groups', '', 'Quality', '');
        CondGroups.Ngroups = length(cfg.channels);
        CondGroups.Groups = num2cell(1:nchannels);
    else
        fprintf('Defining conditioning sets...')
        CondGroups = struct('Permutation', num2cell([1:cfg.nperms]'), 'Ngroups', '', 'Groups', '', 'Quality', '');
        pm = 1;
        while pm <= cfg.nperms
            
            ngroups = randi([cfg.ngroups(1), cfg.ngroups(2)], 1); % choose number of groups)
            CondGroups(pm).Ngroups = ngroups;
            
            nChanGroup = zeros(ngroups,1);
            for nc = 1:ngroups-1 % choose number of channels per group
                nChanGroup(nc) = randi(nchannels-sum(nChanGroup)-ngroups+nc,1);
            end
            nChanGroup(end) = nchannels-sum(nChanGroup);
            
            ChanIDPerm = randperm(nchannels,nchannels);
            
            % make groups
            tmpcs = [0;cumsum(nChanGroup)];
            condSelID = cell(ngroups,1);
            for ng = 1:ngroups
                condSelID{ng} = ChanIDPerm(tmpcs(ng)+1:tmpcs(ng+1));
            end
            
            % check whether cond. set has occured before
            breakloop = 0;
            for perm = 1:pm-1
                %                 cellfun(@isequal, condSelID, CondGroups(perm).Groups)
                if length(condSelID) == length(CondGroups(perm).Groups)
                    if cellfun(@isequal, condSelID, CondGroups(perm).Groups)
                        breakloop = 1;
                    end
                end
            end
            if breakloop == 1
                continue
            else
                CondGroups(pm).Groups = condSelID;
                pm = pm +1;
            end
        end
    end
    fprintf(' Done!\n')
    
    % save conditioning sets for soft restart
    save([cfg.path, 'PartialCoherence/', cfg.subject, ...
        cfg.fileName, '_CondGroups','.mat'], 'CondGroups');
end
end