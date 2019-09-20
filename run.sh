kDim=400
numberEventsToSavePerProcess=30
nProcess=30
seedShift=12314
nentries=500
override_nentries=0
verbose=0

g++ -o main main.C `root-config --cflags --glibs`
rm histograms/*
rm logs/*
rm diagnostic_logs.txt

for ((iProcess=0; iProcess < $nProcess; iProcess++));
do
    ./main "$iProcess" "$kDim" "$numberEventsToSavePerProcess" "$nProcess" "$seedShift" "$nentries" "$override_nentries" "$verbose" &
done

wait

cat logs/* > diagnostic_logs.txt
root -l -b -q makeDiagnosticHists.C
