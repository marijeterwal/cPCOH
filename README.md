# Consensus-based Partial Coherence

Partial coherence was introduced in Rosenberg et al., Journal of Neuroscience Methods, 1998. Like normal coherence, partial coherence captures the frequency-resolved correlation of a pair of interest. In computing the partial coherence, the influences of all other recording sites in the data set are conditioned out, hence removing spurious coherence introduced by shared inputs to the pair of interest, volume conduction, etc. Partial coherence can not be reliably computed for datasets that contain a high number of recording sites or few trials. Instead, conditioning can be performed using averages of randomly grouped recording sites. When this grouping is repeated multiple times a consensus partial coherence can be established, i.e. consensus-based partial coherence.

This Matlab code implements the consensus-based partial coherence method presented in Ter Wal et al., NeuroImage, 2018 (https://doi.org/10.1016/j.neuroimage.2018.06.011)
The code computes the time- and scale-resolved wavelet coherence across trials on FieldTrip-style data structures. It contains a sample dataset and an example script that demonstrates how the code can be used.


## Prerequisites
- Matlab (The MathWorks) - The code has been tested on several Matlab versions between R2014a and R2018a
- Wavelet Toolbox for Matlab (The Mathworks)


## Getting started
1. Download the .zip and unzip or clone to your favourite path.
2. Make sure you Matlab path is set to include the code (Set path on the Home tab).
3. Walk throught the example script. It demonstrates how to set up the config and call the computeCPCOH function.
4. Try it on your own data (NOTE: your data will have to organised as a FieldTrip data struct).


## License and disclaimer
This code is free to use and modify for everyone, but comes without warranty. 


## Citation
If you find this code useful, please refer to the Github repository and cite the paper:
- https://github.com/marijeterwal/cPCOH/
- Ter Wal et al., NeuroImage, 2018. https://doi.org/10.1016/j.neuroimage.2018.06.011

