echo "Deleting old log files"
rm -f log_bcal.txt
rm -f log_fcal.txt
rm -f log_split.txt

echo "Cleaning directories" 
./cleanDirectories.sh

echo "Starting fcal q-value analysis..."
sed -i 's@^detector=".*"@detector="fcal"@g' run.py
python run.py > log_fcal.txt
echo "Starting bcal q-value analysis..."
sed -i 's@^detector=".*"@detector="bcal"@g' run.py
python run.py > log_bcal.txt
echo "Starting split q-value analysis..."
sed -i 's@^detector=".*"@detector="split"@g' run.py
python run.py > log_split.txt

echo "hadding bcal/fcal/split histograms"
rm -f postQVal.root
hadd diagnosticPlots/postQVal.root diagnosticPlots/bcal/postQValHists_bcal.root diagnosticPlots/fcal/postQValHists_fcal.root diagnosticPlots/split/postQValHists_split.root
echo "hadding bcal/fcal/split flatTrees"
rm -f postQVal_flatTree.root
hadd diagnosticPlots/postQVal_flatTree.root diagnosticPlots/bcal/postQ_bcal_flatTree.root diagnosticPlots/fcal/postQ_fcal_flatTree.root diagnosticPlots/split/postQ_split_flatTree.root

echo "Drawing the stacked histograms of the summed bcal/fcal/split datasets"
root -l -b -q makeDiagnosticHists_drawSum.C
