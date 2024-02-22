================================================================================================================================================================================================

CONTENTS OF SCORES
------------------

+ 'userOpinionScores'
      This folder contains the data collected from the participants. The
      data is represented in CSV files, including the name of the stimuli,
      participant vote, presentation time (i.e., time spent while the
      stimulus is shown), and decision time (i.e., time spent after end 
      of presentation until the vote).

- 'analyzeResults.m' 
      This m-file is used to treat the individual user opinion scores as
      described in the paper. When this script is run in Matlab, it 
      computes the mean opinion scores (MOS) and confidence intervals (CI)
      for each stimulus, plots the figures for each content, and saves
      the scores in a mat file (commented out to avoid overwriting).

- 'vsenseVVDB2_bitrates.mat' 
      This file contains the bitrates for all distorted (i.e., compressed) 
      stimuli. It includes following variables:
       == namesVecDist .......: names of the compressed stimuli.
       == fileSizeVecDist ....: file sizes for each stimulus in bits.
       == rateKbpsVecDist ....: bitrate in Kbps for each stimulus.
       == rateMbpsVecDist ....: bitrate in Mbps for each stimulus.

- 'vsenseVVDB2_subjectiveScores.mat' 
      This file contains the subjective quality scores obtained. It 
      includes following variables:
       == indsOrig ...........: indices for reference stimuli.
       == resultStr ..........: a struct for all the results:
            -- wRef_names     : stimuli names including references.
            -- wRef_MOS       : MOS for each stimuli. (inc. refs.)
            -- wRef_CI        : CI for each stimuli. (inc. refs.)
            -- Dist_names     : stimuli names excluding references.
            -- Dist_compTypes : compression types for compressed stimuli.
            -- Dist_compTypeQ : list of compression types. (4x1 vector)
            -- Dist_contentQ  : list of vol. vid. contents. (8x1 vector)
            -- Dist_MOS       : MOS for each stimuli. (only dist.)
            -- Dist_CI        : CI for each stimuli. (only dist.)

- 'vsenseVVDB2_subjectiveScores.xlsx' 
      This file also contains the subjective quality scores obtained. It 
      includes two sheets: one for all stimuli including references and
      another for all the compressed stimuli (excluding references).

================================================================================================
================================================================================================