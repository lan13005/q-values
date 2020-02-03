rm -r $grid13/graphs/diagnosticPlots
rm -r $grid13/graphs/fitResults
rm -r $grid13/graphs/histograms
rsync -r diagnosticPlots $grid13/graphs
rsync -r fitResults $grid13/graphs
rsync -r histograms $grid13/graphs

