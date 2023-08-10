#! /bin/bash
# TrackEff094
# diffSys = 'TrackEff094', 'BKGNoFit', 'BKGV2', 'V0C', 'V0A'
# source runHaddScaleFiles.sh 'LHC18r' 5  diffSys

for i in `seq 0 4`; do
    echo $i $1 $2 $3
    python MergeScalePtHardPlots.py $i $1 $2 $3
done

python HAddOutputScaledFiles.py $1 $2 $3

# rm -rf ~/ALICE/cernbox/SWAN_projects/outputFiles/$1/pass3/Ch/Embedding/EmbedPtHardScaledResults_TrackPtCut5_*_Cent*

echo root ~/ALICE/cernbox/SWAN_projects/outputFiles/$1/pass3/Ch/Embedding/


# ~/ALICE/cernbox/SWAN_projects/outputFiles/LHC18qr/pass3/Ch/Embedding

# hadd ~/ALICE/cernbox/SWAN_projects/outputFiles/LHC18qr/pass3/Ch/Embedding/EmbedPtHardScaledResults_TrackPtCut5__Ver1.root ~/ALICE/cernbox/SWAN_projects/outputFiles/LHC18q/pass3/Ch/Embedding/EmbedPtHardScaledResults_TrackPtCut5__Ver1.root ~/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Embedding/EmbedPtHardScaledResults_TrackPtCut5__Ver1.root

# hadd ~/ALICE/cernbox/SWAN_projects/outputFiles/LHC18qr/pass3/Ch/Embedding/EmbedPtHardScaledResults_TrackPtCut5_TrackEff094_Ver1.root ~/ALICE/cernbox/SWAN_projects/outputFiles/LHC18q/pass3/Ch/Embedding/EmbedPtHardScaledResults_TrackPtCut5_TrackEff094_Ver1.root ~/ALICE/cernbox/SWAN_projects/outputFiles/LHC18r/pass3/Ch/Embedding/EmbedPtHardScaledResults_TrackPtCut5_TrackEff094_Ver1.root


