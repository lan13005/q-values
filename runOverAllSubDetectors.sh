echo "Deleting old log files"
rm -f log_bcal.txt
rm -f log_fcal.txt
rm -f log_split.txt

echo "Starting fcal q-value analysis..."
sed -i 's@detector=".*"@detector="fcal"@g' main.C
sed -i 's@detector=".*"@detector="fcal"@g' run.py
sed -i 's@detector=".*"@detector="fcal"@g' makeDiagnosticHists.C
python run.py > log_fcal.txt
echo "Starting bcal q-value analysis..."
sed -i 's@detector=".*"@detector="bcal"@g' main.C
sed -i 's@detector=".*"@detector="bcal"@g' run.py
sed -i 's@detector=".*"@detector="bcal"@g' makeDiagnosticHists.C
python run.py > log_bcal.txt
echo "Starting split q-value analysis..."
sed -i 's@detector=".*"@detector="split"@g' main.C
sed -i 's@detector=".*"@detector="split"@g' run.py
sed -i 's@detector=".*"@detector="split"@g' makeDiagnosticHists.C
python run.py > log_split.txt


echo "hadding bcal/fcal/split histograms"
rm -f postQVal.root
hadd postQVal.root postQValHists_bcal.root postQValHists_fcal.root postQValHists_split.root
echo "hadding bcal/fcal/split flatTrees"
rm -f postQVal_flatTree.root
hadd postQVal_flatTree.root postQ_bcal_flatTree.root postQ_fcal_flatTree.root postQ_split_flatTree.root

echo "Drawing the stacked histograms of the summed bcal/fcal/split datasets"
root -l -b -q makeDiagnosticHists_drawSum.C
