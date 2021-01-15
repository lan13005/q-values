for i in {0..10..1}
do
    root -l -b -q "checkGausFeedThrough.C($i)" > checkGausFeedThrough/log2_$i.txt &
done
