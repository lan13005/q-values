start_time="$(date -u +%s)"

kDim=200
numberEventsToSavePerProcess=1
nProcess=30
seedShift=125
nentries=100000
override_nentries=1
verbose=0
varString="cosTheta_X_cms;phi_X_cms;cosTheta_eta_gjs;phi_eta_gjs" #vanHove_omegas"


numVars=$(($(grep varString run.sh | grep -o ";" | wc -l)))
sed -i "s@dim=dimNum@dim=$numVars@g" main.h

g++ -o main main.C `root-config --cflags --glibs`

sed -i "s@dim=$numVars@dim=dimNum@g" main.h
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


end_time="$(date -u +%s)"
elapsed="$(($end_time-$start_time))"
echo "Total of $elapsed seconds elapsed for process"
