% create CPCOHstruct

% Description:
% This function is generates a data struct to store COH and cPCOH results.

% INPUT: 
% The code requires four input arguments:
    % cfg: a config struct specifying the main parameters that have been set
    % data: a FieldTrip data struct with preprocessed data;
    % pair: array with IDs of the channels in the pair
    % perm: scalar or array of permutation IDs

% OUTPUT: 
% The functions outputs the emtpy data struct outStruct.

% This code implements the consensus-based partial coherence method as
% described in Ter Wal et al., NeuroImage, 2018.
% DOI: https://doi.org/10.1016/j.neuroimage.2018.06.011
% For more details and citations please refer to the paper. 

% Marije ter Wal, 2018
% m.terwal@donders.ru.nl || m.j.terwal@bham.ac.uk

function outStruct = createCPCOHstruct(cfg, data, pair, perm)

nfreqs      = length(cfg.foi);
ntimes      = length(cfg.toi);
nevents     = length(cfg.events);

outStruct = rmfield(data,{'trial','time'});
outStruct.time = cfg.toi;
outStruct.freq = cfg.foi;
outStruct.ev = cfg.events;
outStruct.perm = perm;
outStruct.dimord = 'event_perm_freq_time';
outStruct.cfg = cfg;
outStruct.label ={data.label{pair(1)},data.label{pair(2)}};
[outStruct.Coh, outStruct.CohZ, outStruct.PCoh, outStruct.PCohZ] = ...
    deal(zeros(nevents,length(perm),nfreqs,ntimes));
if length(perm) > 1
    [outStruct.Consensus, outStruct.PCohc, outStruct.PCohZc] = ...
    deal(zeros(nevents,1,nfreqs,ntimes));
end

end