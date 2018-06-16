% checkConfig

% Description:
% This function is called by getSpectra and copmuteCPCOH to check and set 
% the main parameters of the consensus-based partial coherence.

% INPUT: 
% The code requires three input arguments:
    % cfg: a config struct specifying the main parameters that have been set
    % data: a FieldTrip data struct with preprocessed data;
    % type: string to indicate what settings have to be checeked ('SPECTRA', 'PCOH').

% OUTPUT: 
% The functions outputs the checked cfg struct.

% This code implements the consensus-based partial coherence method as
% described in Ter Wal et al., NeuroImage, 2018.
% DOI: https://doi.org/10.1016/j.neuroimage.2018.06.011
% For more details and citations please refer to the paper. 

% Marije ter Wal, 2018
% m.terwal@donders.ru.nl || m.j.terwal@bham.ac.uk


function cfg = checkConfig(cfg, data, type)

if ~isfield(cfg, 'path');               cfg.path = data.cfg.path; end
if ~isfield(cfg, 'subject');            cfg.subject = data.cfg.subject; end

switch type
    
    case 'SPECTRA'
        % general
        if ~isfield(cfg, 'channels');       cfg.channels = 1:length(data.label); end
        if ~isfield(cfg, 'Zscore');         cfg.Zscore = false; end

        % frequencies
        if ~isfield(cfg,'df');              cfg.df = 1; end
        if ~isfield(cfg,'foi') && ~isfield(cfg,'foilim')
                                            error('No frequencies specified')
        elseif isfield(cfg, 'foilim')
                                            cfg.foi = cfg.foilim(1):cfg.df:cfg.foilim(2);
        end
        cfg.scales                          = data.fsample./cfg.foi;
        if ~isfield(cfg,'wavelet');         cfg.wavelet = 'cmor3-1'; end
        
      % times
        if ~isfield(cfg, 'dt')          
            if isfield(cfg, 'toi') && length(cfg.toi) > 2
                                            cfg.dt = cfg.toi(2)-cfg.toi(1);
            else
                                            cfg.dt = 1/data.fsample; 
            end
        end
        if ~isfield(cfg,'toi') && ~isfield(cfg,'toilim')
                                            error('No times specified')
        elseif isfield(cfg, 'toilim')
                                            cfg.toi = cfg.toilim(1):cfg.dt:cfg.toilim(2);
        end    
        
    case 'PCOH'
        % general
        if ~isfield(cfg,'channels');        cfg.channels = 1:length(data.label); end
        if ~isfield(cfg,'pairs');           cfg.pairs = combnk(cfg.channels,2); end
        if ~isfield(cfg,'events');          cfg.events = 1:size(data.event,1); end
        
        if ~isfield(cfg,'overwrite');       cfg.overwrite = false; end
        if ~isfield(cfg,'deleteTemp');      cfg.deleteTemp = false; end
        if ~isfield(cfg,'alpha');           cfg.alpha = 0.05; end
        if ~isfield(cfg,'consensusThres');  cfg.consensusThres = 0.9; end
        if ~isfield(cfg,'weightFunction');  cfg.weightFunction = @(x) heaviside(x-cfg.consensusThres); end
        
        % frequencies
        if ~isfield(cfg,'df');              cfg.df = 1; end
        if ~isfield(cfg,'foi') && ~isfield(cfg,'foilim')
            error('No frequencies specified')
        elseif isfield(cfg, 'foilim')
            cfg.foi = cfg.foilim(1):cfg.df:cfg.foilim(2);
        end
        cfg.scales                          = data.fsample./cfg.foi;
        if ~isfield(cfg,'wavelet');         cfg.wavelet = 'cmor3-1'; end
        
        % times
        if ~isfield(cfg, 'dt')          
            if isfield(cfg, 'toi') && length(cfg.toi) > 2
                cfg.dt = cfg.toi(2)-cfg.toi(1);
            else
                cfg.dt = 1/data.fsample; 
            end
        end
        if ~isfield(cfg,'toi') && ~isfield(cfg,'toilim')
                                            error('No times specified')
        elseif isfield(cfg, 'toilim')
                                            cfg.toi = cfg.toilim(1):cfg.dt:cfg.toilim(2);
        end
        
        % partial coherence - general
        if ~isfield(cfg,'nperms');          cfg.nperms = 50; end
        if ~isfield(cfg,'ngroups');         cfg.ngroups = [2,4]; end
        
        % determine type of reference
        if isfield(cfg, 'baseline') % baseline in the same dataset
            cfg.nreps = length(cfg.baseline(1):cfg.dt:cfg.baseline(2));
        elseif ~isfield(cfg, 'nreps')
            cfg.nreps = size(dataRef.time{1});
        end
        
end
end