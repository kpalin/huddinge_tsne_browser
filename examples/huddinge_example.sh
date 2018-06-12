#!/bin/bash

#  Tue Jun 12 09:00:27 2018
#  kpalin@merit-ltdk.it.helsinki.fi

# For debugging
#set -o verbose 

# Die on unset variables
set -o nounset
# Die on errors
set -o errexit
# Die if any part of a pipe fails
set -o pipefail


echo "Download distances"

wget --no-clobber http://www.cs.helsinki.fi/u/kpalin/all8mers_min_rev_complement.dists.gz
gunzip all8mers_min_rev_complement.dists.gz || true

echo "Download selex data"
wget --no-clobber ftp.sra.ebi.ac.uk/vol1/fastq/ERR100/006/ERR1003746/ERR1003746.fastq.gz
wget --no-clobber ftp.sra.ebi.ac.uk/vol1/fastq/ERR100/008/ERR1003748/ERR1003748.fastq.gz
wget --no-clobber ftp.sra.ebi.ac.uk/vol1/fastq/ERR100/000/ERR1003750/ERR1003750.fastq.gz
wget --no-clobber ftp.sra.ebi.ac.uk/vol1/fastq/ERR100/002/ERR1003752/ERR1003752.fastq.gz
wget --no-clobber ftp.sra.ebi.ac.uk/vol1/fastq/ERR100/004/ERR1003874/ERR1003874.fastq.gz
wget --no-clobber ftp.sra.ebi.ac.uk/vol1/fastq/ERR100/006/ERR1003876/ERR1003876.fastq.gz
wget --no-clobber ftp.sra.ebi.ac.uk/vol1/fastq/ERR100/008/ERR1003878/ERR1003878.fastq.gz
wget --no-clobber ftp.sra.ebi.ac.uk/vol1/fastq/ERR100/000/ERR1003880/ERR1003880.fastq.gz

K=8
echo "Calculating ${K}mers"

for i in *.fastq.gz
do
    
    OUT=$(basename $i .fastq.gz).${K}mer_counts.jf
    if [ ! -e ${OUT} ];
    then
        zcat $i | jellyfish count --canonical -o $OUT --text -m ${K} -s 1M --bf-size 1G -t 16 --disk /dev/stdin &
    fi
done
wait


echo "Make config files"
cat >HNF4_config_example.json <<EOF
{
   "_config": {
      "kmer_size": 8
   }, 
   "HNF4A": [
      {
         "fasta": "ERR1003746.fastq.gz", 
         "filename": "ERR1003746.8mer_counts.jf", 
         "cycle": 1
      }, 
      {
         "fasta": "ERR1003748.fastq.gz", 
         "filename": "ERR1003748.8mer_counts.jf", 
         "cycle": 2
      }, 
      {
         "fasta": "ERR1003750.fastq.gz", 
         "filename": "ERR1003750.8mer_counts.jf", 
         "cycle": 3
      }, 
      {
         "fasta": "ERR1003752.fastq.gz", 
         "filename": "ERR1003752.8mer_counts.jf", 
         "cycle": 4
      }
   ]
}
EOF

cat >HOXB13_config_example.json <<EOF
{
   "_config": {
      "kmer_size": 8
   }, 
   "HOXB13": [
      {
         "fasta": "ERR1003874.fastq.gz", 
         "filename": "ERR1003874.8mer_counts.jf", 
         "cycle": 1
      }, 
      {
         "fasta": "ERR1003876.fastq.gz", 
         "filename": "ERR1003876.8mer_counts.jf", 
         "cycle": 2
      }, 
      {
         "fasta": "ERR1003878.fastq.gz", 
         "filename": "ERR1003878.8mer_counts.jf", 
         "cycle": 3
      }, 
      {
         "fasta": "ERR1003880.fastq.gz", 
         "filename": "ERR1003880.8mer_counts.jf", 
         "cycle": 4
      }
   ]
}

EOF


echo "Launching image servers"
huddinge_tsne_browser -i all8mers_min_rev_complement.dists -j HNF4_config_example.json  -V &
huddinge_tsne_browser -i all8mers_min_rev_complement.dists -j HOXB13_config_example.json  -V &

echo "Waiting"
wait
