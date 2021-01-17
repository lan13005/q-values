for i in {0..10..1}
do
    root -l -b -q "checkGausFeedThrough.C($i)" > newLogs/log_$i.txt &
done
