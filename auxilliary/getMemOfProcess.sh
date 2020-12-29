id=$1
pid=$2
maxIters=$3
sleepTime=$4
outDir=$5

echo "memory" > $outDir/times_$1.log

for i in $(seq 0 $maxIters)
do
    sleep $sleepTime
    # check if the process is still running. If it is we will continue to output
    if ! kill -s 0 $pid > /dev/null 2>&1
    then
        break
    fi
    ps -o rss $pid | tail -n1 >> $outDir/times_$1.log
done
