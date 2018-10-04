./bin/gtf
sleep 1
./bin/ffa
sleep 1
echo "please ignore the above output. Comparing the from ./output_graphtv.txt and ./output_ffa.txt"
diff output_graphtv.txt output_ffa.txt
