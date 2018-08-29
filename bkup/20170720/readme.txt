——————————————————————
Meta scripts:
Individuals- analysisController.m
Group- groupAnalysisController.m

Parameter files:
Individuals- analysisParams
Group- groupAnalysisParams
——————————————————————


This analysis pipeline processes SSVEPs in EEG data measured during a binocular rivalry task. It performs the following analysis functions:

-FFT analysis
-RLS analysis

Individual Analysis:
-FFT analysis of oscillation rate
-Transition analysis
-Peak analysis
-Participant plots of transitions and peaks
-Amplitude cluster analysis
Functions:


Group Analysis:
- Correlations with behavioral indices
- Group transition and peak plotting
- Group amplitude cluster analysis
- Group comparisons of FFT oscillation analysis
- Test-retest for FFT oscillation analysis
Functions: avgFFT.m, FFTconsistency.m

Files needed prior to analysis:
tLists_CERW/parName_runName.mat (switch data from Caroline)