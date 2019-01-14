The actual data comparison results are based on BWA-mem. In order to process the bam file, 
we first use Samtools to convert bam to sam format. The specific commands are as follows:
    samtools view *.bam > *.sam
After getting the sam file, we use "change_sam" to process it into SVseq-ft format, the specific 
command is as follows:
    python read_sam.py -i <sam_file> -o <output_file>

To facilitate real data test parameters, we divide the SVseq-ft algorithm into two parts:
1.get_score
  python connect.py -i <sam_file> - w <win_size> -o <out_name>
  Through this step we can get the variation score file corresponding to the sam file.
2.filter
  python connect.py -t <thre_value> - v <var_value> -m <min_support> - i <sam_file> - s <score_file> -o <out_file>
  Through this step, we can get the location information and length information of the detected complex indel.

We performed experiments based on sequencing data from NA12878 individual, 
which is published at ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/integrated_sv_map/supporting/NA12878/.