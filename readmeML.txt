%%%
Meta scripts:
groupAnalysisController.m

Parameter files (Required for all scripts below):
svmParams.m
groupAnalysisParams.m
%%%

%%-------------------------------------------------------------------------

%%Set up scripts (uses left out people)
runSVMcrossVal.m - 
runSVMAirplaneMatrices.m - 

%%-------------------------------------------------------------------------

%%Running scripts
runSVMclassifyPercept.m
classSVMAvsC.m -

%%-------------------------------------------------------------------------

%%Classifying scripts
SVMclassify.m
SVMclassify_withMixed.m 
-- classifies 
-- called by runSVMclassifyPercept.m
SVMclassify_leftOutPar.m

%%-------------------------------------------------------------------------

%%Scripts that build training and testing sets

makeSVMmatrices.m 
-- called by: svmClassify.m, svmClassify_mixed.m, svmClassLeftOutPar.m, runSVMcrossVal.m
-- dependencies: getSVMSegs.m, getSVMSegPresses.m
-- Creates matrices of predictors (each column is a feature, each row is an 
   example) [trainingData, testData] and labels [trainingLabels, testLabels],
   to be used as training and testing sets
 
getSVMSegs.m 
-- called by: makeSVMmatrices.m
-- dependencies: none
-- Collects matrices of segments (each row is a segment) for the low 
  [lowSegList] and high [highSegList] frequency traces. labelVect is a cell
   array with percept labels (hi lo and mi) corresponding to each segment

getSVMTmatrix.m 
-- called by: 
-- dependencies:

getSVMSegPresses.m	- called by XXX,

%%-------------------------------------------------------------------------

%%Plotting scripts
plotSVMtimeCourses.m - 
plotSVMAirplaneMat.m - 

%%-------------------------------------------------------------------------

