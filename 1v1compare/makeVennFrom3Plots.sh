cat $1 | ./extract.sh > out1
cat $2 | ./extract.sh > out2
cat $3 | ./extract.sh > out3

Rscript ./makeVenns.R out1 out2 out3
