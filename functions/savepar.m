% savepar

% Description:
% This function makes the save function compatible with the parfor loop.

% This code belongs to the consensus-based partial coherence method as
% described in Ter Wal et al., NeuroImage, 2018.
% DOI: https://doi.org/10.1016/j.neuroimage.2018.06.011
% For more details and citations please refer to the paper. 

% Marije ter Wal, 2018
% m.terwal@donders.ru.nl || m.j.terwal@bham.ac.uk


function savepar(path, var, vname)

eval([vname,'= var;']);
save(path, vname);

end