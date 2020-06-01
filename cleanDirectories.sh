rm main

rm -r logs
mkdir -p logs/bcal
mkdir -p logs/fcal
mkdir -p logs/split

rm -r histograms
mkdir -p histograms/fcal
mkdir -p histograms/bcal
mkdir -p histograms/split

rm -r diagnosticPlots
mkdir -p diagnosticPlots/fcal
mkdir -p diagnosticPlots/bcal
mkdir -p diagnosticPlots/split

rm -r fitResults
mkdir fitResults

rm -f log_bcal.txt
rm -f log_fcal.txt
rm -f log_split.txt
