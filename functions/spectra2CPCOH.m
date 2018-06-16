% spectra2CPCOH

% Description:
% This function computes the auto- and cross-spectra, computes group averages
% and performs conditioning of coherence using those averages. 

% INPUT: 
% The code requires two input arguments:
    % data: time-resolved wavelet coefficients for one scale as produces by
    % getSpectra;
    % pairSel: IDs of the channels in the pair of interest;
    % condSel: conditioning group details for one permutation as produces
    % by getPermutations.

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


function [coh, pcoh, qtok] = spectra2CPCOH(data, pairSel, condSel)

qtok = 0;

ngroups     = length(condSel); % number of groups to condition on
ntimes      = size(data,3);

% preallocate
pcoh = nan(ntimes,1,'single');
coh = nan(ntimes,1,'single'); %0;

% determine partial coherence for each pair
nc1 = pairSel(1);
nc2 = pairSel(2);

%% Prepare conditioning dataset by averaging within groups

dataCond = cat(2,data(:,nc1,:), data(:,nc2,:));
for ng = 1:ngroups
    pch = setdiff(condSel{ng}, [nc1,nc2]); % channels to condition on
    if isempty(pch); ngroups = ngroups-1; continue; end
    tmp = mean(data(:,pch,:),2);
    dataCond = cat(2,dataCond,tmp);
end

if size(dataCond,2) == 2
    fprintf('\nNot a valid conditioning set, conditioning aborted.')
    return
end

%% Compute averaged cross-spectra matrix 
% This has to be a full and symmetric matrix of all auto- and
% cross-spectra.

ACond = zeros(ngroups+2,ngroups+2, ntimes, 'single');

for cch1 = 1:ngroups+2
    for cch2 = cch1:ngroups+2
		% compute the upper triangle of the matrix
        ACond(cch1,cch2,:) = mean(abs( dataCond(:,cch1,:).*conj(dataCond(:,cch2,:)) ),1);
        % complete the lower triangle of the matrix
        if cch1 ~= cch2
            ACond(cch2,cch1,:) = conj(ACond(cch1,cch2,:));
        end
    end
end

%% Conditioning for each point in time 

for t = 1:ntimes
    % unconditioned auto- and cross-spectra 
    A11 = ACond(1,1,t);
    A22 = ACond(2,2,t);
    A12 = ACond(1,2,t);
    
    % inverse of conditioning matrix
    BB = inv(ACond(3:end,3:end,t)); 
	% see for more details on the LU decomposition and SVD based inverses:
	% Press, W. H., Teukosly, S. A., Vetterling, W. T., & Flannery, B. P. (2007). 
	% Numerical Recipes. The Art of Scientific Computing. (3rd Ed). Cambridge Univ Press, Cambridge, UK
	
    % conditioned auto- and cross-spectra 
    Cp12 = A12 - ACond(1,3:end,t) * BB * ACond(3:end,2,t);
    Cp11 = A11 - ACond(1,3:end,t) * BB * ACond(3:end,1,t);
    Cp22 = A22 - ACond(2,3:end,t) * BB * ACond(3:end,2,t);
    
    % compute COH and cPCOH
    coh(t) = abs(A12) / sqrt(A11 * A22); % wavelet coherence
    pcoh(t) = abs(Cp12) / sqrt(Cp11 * Cp22); % partial wavelet coherence
end

qtok = 1;

end






