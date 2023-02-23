# Neuronal-Network-Inference-with-Reservoir-Computer
Code to build and train a reservoir computer to infer neuronal connectivity from calcium fluorescence time series data of *C. elegans* neurons.

This repository contains a MATLAB code (elegans_network_inference.m) demonstrating a minimal example of a Reservoir Computer link scoring technique. The code is associated with the paper "Network inference from short, noisy, low time-resolution, partial measurements: Application to *C. elegans* neuronal calcium dynamics" by Amitava Banerjee, Sarthak Chandra, and Edward Ott, published in the Proceedings of the National Academy of Sciences (PNAS).

Source of the time series dataset: https://osf.io/na4f9, Whole brain imaging data from Kato et al., Cell (2015): Global Brain Dynamics Embed the Motor Command Sequence of *Caenorhabditis elegans*, https://pubmed.ncbi.nlm.nih.gov/26478179/. **This repository does not contain the time series dataset. The time series dataset must be downloaded from the database at https://osf.io/na4f9 separately by the user, and must be stored with the file name WT_NoStim.mat (the original file name) in the same directory as the code, for the code to run.**

Source of the ground truth neuronal connectivity of *C. elegans*: https://www.wormatlas.org/neuronalwiring.html#Connectivitydata . We folded over the left-right neuron pairs as described in the paper. The connectivity is stored as a binary matrix J0.mat in the same directory as the code. The matrix is 8-by-8, corresponding to 8 neurons in the order ["AVA","SMDV","RIV","RIM","SMDD","AIB","RIB","AVB"].

The main MATLAB code is elegans_network_inference.m . Please put the code in the same directory as the datafiles WT_NoStim.mat and J0.mat . The code, upon running, plots (1) the time-series data for the above 8 pairs of neurons, (2) the actual and inferred connectivity matrices between the neurons, and (3) histograms of reservoir computer inferred scores for neurons which are connected and neurons which are not.

Contact: Amitava Banerjee at amitava8196@gmail.com or amitavab@cshl.edu for questions about the code.
