kDim=200
numberEventsToSavePerProcess=10
nProcess=40
seedShift=12314
nentries=300
override_nentries=0
verbose=0
varString="cosTheta_X_cms;phi_X_cms;cosTheta_eta_gjs;phi_eta_gjs;cosThetaHighestEphotonIneta_gjs;cosThetaHighestEphotonInpi0_cms"

g++ -o main main.C `root-config --cflags --glibs`
rm histograms/*
rm logs/*
rm diagnostic_logs.txt

for ((iProcess=0; iProcess < $nProcess; iProcess++));
do
    ./main "$iProcess" "$kDim" "$numberEventsToSavePerProcess" "$nProcess" "$seedShift" "$nentries" "$override_nentries" "$verbose" $varString &
done

wait

cat logs/log* > diagnostic_logs.txt
rm qvalResults.root
hadd qvalResults.root logs/results*
root -l -b -q makeDiagnosticHists.C
