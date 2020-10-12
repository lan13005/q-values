startTime=$SECONDS

condor_submit submit.main
# since bash's for loop is inclusive we have to subtract 1 from the total counts
counts=$((36-1))
for i in $(eval echo {0..$counts} )
do
    condor_wait condor/job$i/test.log
done

endTime=$SECONDS
echo "elapsed time: $((endTime-startTime)) seconds"
