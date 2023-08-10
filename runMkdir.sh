#! /bin/bash

for iCent in `seq 0 4`; do
    mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18q/pass3/Ch/Unfolded/UnfoldPerformancePlot//InPlane/Cent$iCent/Bayes/RatioPlot
    mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Unfolded/UnfoldPerformancePlot//InPlane/Cent$iCent/Bayes/RatioPlot
    mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18qr/pass3/Ch/Unfolded/UnfoldPerformancePlot//InPlane/Cent$iCent/Bayes/RatioPlot
    mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18q/pass3/Ch/Unfolded/UnfoldPerformancePlot//OutOfPlane/Cent$iCent/Bayes/RatioPlot
    mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Unfolded/UnfoldPerformancePlot//OutOfPlane/Cent$iCent/Bayes/RatioPlot
    mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18qr/pass3/Ch/Unfolded/UnfoldPerformancePlot//OutOfPlane/Cent$iCent/Bayes/RatioPlot
    for iIte in `seq 1 10`; do
        echo $iIte
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18q/pass3/Ch/Unfolded/UnfoldPerformancePlot//InPlane/Cent$iCent/Bayes/Iteration$iIte/CorrelationCoefficients/
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18q/pass3/Ch/Unfolded/UnfoldPerformancePlot//InPlane/Cent$iCent/Bayes/Iteration$iIte/FoldedPbPbTruth
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18q/pass3/Ch/Unfolded/UnfoldPerformancePlot//InPlane/Cent$iCent/Bayes/Iteration$iIte/RefoldingTest
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18q/pass3/Ch/Unfolded/UnfoldPerformancePlot//InPlane/Cent$iCent/Bayes/Iteration$iIte/ClosureTest
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18q/pass3/Ch/Unfolded/UnfoldPerformancePlot//InPlane/Cent$iCent/Bayes/Iteration$iIte/UnfoldedMeasuredRatio
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18q/pass3/Ch/Unfolded/UnfoldPerformancePlot//InPlane/Cent$iCent/Bayes/Iteration$iIte/KinEff

        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18q/pass3/Ch/Unfolded/UnfoldPerformancePlot//OutOfPlane/Cent$iCent/Bayes/Iteration$iIte/CorrelationCoefficients/
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18q/pass3/Ch/Unfolded/UnfoldPerformancePlot/OutOfPlane/Cent$iCent/Bayes/Iteration$iIte/FoldedPbPbTruth
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18q/pass3/Ch/Unfolded/UnfoldPerformancePlot/OutOfPlane/Cent$iCent/Bayes/Iteration$iIte/RefoldingTest
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18q/pass3/Ch/Unfolded/UnfoldPerformancePlot/OutOfPlane/Cent$iCent/Bayes/Iteration$iIte/ClosureTest
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18q/pass3/Ch/Unfolded/UnfoldPerformancePlot/OutOfPlane/Cent$iCent/Bayes/Iteration$iIte/UnfoldedMeasuredRatio
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18q/pass3/Ch/Unfolded/UnfoldPerformancePlot/OutOfPlane/Cent$iCent/Bayes/Iteration$iIte/KinEff

        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Unfolded/UnfoldPerformancePlot//InPlane/Cent$iCent/Bayes/Iteration$iIte/CorrelationCoefficients/
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Unfolded/UnfoldPerformancePlot//InPlane/Cent$iCent/Bayes/Iteration$iIte/FoldedPbPbTruth
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Unfolded/UnfoldPerformancePlot//InPlane/Cent$iCent/Bayes/Iteration$iIte/RefoldingTest
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Unfolded/UnfoldPerformancePlot//InPlane/Cent$iCent/Bayes/Iteration$iIte/ClosureTest
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Unfolded/UnfoldPerformancePlot//InPlane/Cent$iCent/Bayes/Iteration$iIte/UnfoldedMeasuredRatio
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Unfolded/UnfoldPerformancePlot//InPlane/Cent$iCent/Bayes/Iteration$iIte/KinEff

        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Unfolded/UnfoldPerformancePlot//OutOfPlane/Cent$iCent/Bayes/Iteration$iIte/CorrelationCoefficients/
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Unfolded/UnfoldPerformancePlot/OutOfPlane/Cent$iCent/Bayes/Iteration$iIte/FoldedPbPbTruth
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Unfolded/UnfoldPerformancePlot/OutOfPlane/Cent$iCent/Bayes/Iteration$iIte/RefoldingTest
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Unfolded/UnfoldPerformancePlot/OutOfPlane/Cent$iCent/Bayes/Iteration$iIte/ClosureTest
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Unfolded/UnfoldPerformancePlot/OutOfPlane/Cent$iCent/Bayes/Iteration$iIte/UnfoldedMeasuredRatio
        # mkdir /Users/tkumaoka/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Unfolded/UnfoldPerformancePlot/OutOfPlane/Cent$iCent/Bayes/Iteration$iIte/KinEff
    done
done