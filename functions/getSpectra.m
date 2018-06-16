% getSpectra

% Description:
% This function is called by computeCPCOH to compute the wavelet transform.

% DEPENDENCIES:
% This function uses the cwt function from the Wavelet Toolbox (The
% MathWorks). The syntax from the old (pre-2016b) version is used, but also
% works for the new version of the cwt-function.

% INPUT:
% The code requires two input arguments:
    % cfg: a config struct specifying the main parameters that have been set
    % data: a FieldTrip data struct with preprocessed data;

% OUTPUT:
% The functions outputs the variable cfs: an array of size
% nevents x ntrials x nchannels x nscales x ntimes.

% This code implements the consensus-based partial coherence method as
% described in Ter Wal et al., NeuroImage, 2018.
% DOI: https://doi.org/10.1016/j.neuroimage.2018.06.011
% For more details and citations please refer to the paper.

% Marije ter Wal, 2018
% m.terwal@donders.ru.nl || m.j.terwal@bham.ac.uk


function cfs = getSpectra(cfg, data)

cfg = checkConfig(cfg, data, 'SPECTRA');

% cfs: trials x channels x events x freqs x times
cfs = zeros(length(cfg.events),length(data.trial),length(cfg.channels),length(cfg.foi),length(cfg.toi), 'single');

for tr = 1:length(data.trial)
    for ch = 1:length(cfg.channels)
        fulltrial = cwt(data.trial{tr}(cfg.channels(ch),:), cfg.scales, cfg.wavelet);
        
        % select event and downsample
        for ev = 1:length(cfg.events)
            [~,eventID] = min(abs(data.time{tr}-data.event{ev,tr}));
            TimeIDs = eventID + round(cfg.toi*data.fsample);
            cfs(ev,tr,ch,:,:) = fulltrial(:,TimeIDs);
        end
    end
end

end
