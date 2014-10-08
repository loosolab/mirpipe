mirpipe.pl -file ./data/test.fastq -ref ./data/mirbase20_mature.fa
DIFF=$(diff target.output mirpipe_mirna.tsv) 
if [ "$DIFF" != "" ] 
then
    echo "Automatic test failed. Please check the content of mirpipe_mirna.tsv manually."
fi
if [ "$DIFF" == "" ]
then
    echo "Automatic testing was successful. You are ready to use MIRPIPE."
fi

