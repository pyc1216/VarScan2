# VarScan2
version=v2.3.9.a.1

1. Modify function **comparePositions()** in *Somatic.java*

Use **normalMinVarFreqDecideGenotype** instead of **MinVarFreq** to decide Normal's genotype when Tumor and Normal are both HET(**CASE 2E**).

**normalMinVarFreqDecideGenotype** = 0.03, when **MinVarFreq** < 0.03, else = **MinVarFreq**.

It make sense when **MinVarFreq** is pretty low(e.g 0.005), avoid taking normal HET falsely.

This may call **MORE** somatic snp/indel!


The difference between my version and origin is showed below:
test_data/test.snp.vcf is my version's result, the other is origin version's.

`diff --unified=0 test_data/test.snp.vcf /data_sas/ypu/git_repository/Med1CDx/test_data/ZY/call/chrom/sample.chr7\:10000-159128663.snp.vcf`

    --- test_data/test.snp.vcf	2019-03-08 11:41:03.059672490 +0800
    +++ /data_sas/ypu/git_repository/Med1CDx/test_data/ZY/call/chrom/sample.chr7:10000-159128663.snp.vcf  2019-03-05 16:02:08.653294173 +0800
    @@ -1871 +1871 @@
    -chr7	140453136	.	A	T	.	PASS	DP=2279;SOMATIC;SS=2;SSC=255;GPV=1E0;SPV=8.1581E-59	GT:GQ:DP:RD:AD:FREQ:DP4	0/1:.:378:372:6:1.59%:258,114,5,1	0/1:.:1901:1189:712:37.45%:1001,188,636,76
    +chr7	140453136	.	A	T	.	PASS	DP=2279;SS=1;SSC=255;GPV=1E0;SPV=8.1581E-59	GT:GQ:DP:RD:AD:FREQ:DP4	0/1:.:378:372:6:1.59%:258,114,5,1	0/1:.:1901:1189:712:37.45%:1001,188,636,76
    
    @@ -2672 +2672 @@
    -chr7	151933014	.	C	T	.	PASS	DP=1922;SS=1;SSC=16;GPV=1.7461E-25;SPV=2.5102E-2	GT:GQ:DP:RD:AD:FREQ:DP4	0/1:.:383:372:9:2.36%:209,163,1,8	0/1:.:1539:1466:72:4.68%:1242,224,4,68
    +chr7	151933014	.	C	T	.	PASS	DP=1922;SS=1;SSC=16;GPV=1E0;SPV=2.5102E-2	GT:GQ:DP:RD:AD:FREQ:DP4	0/1:.:383:372:9:2.36%:209,163,1,8	0/1:.:1539:1466:72:4.68%:1242,224,4,68
    
    @@ -2681 +2681 @@
    -chr7	151935689	.	T	C	.	PASS	DP=1586;SS=1;SSC=16;GPV=2.2796E-8;SPV=2.4017E-2	GT:GQ:DP:RD:AD:FREQ:DP4	0/1:.:458:201:4:1.95%:196,5,0,4	0/1:.:1128:346:21:5.72%:342,4,0,21
    +chr7	151935689	.	T	C	.	PASS	DP=1586;SS=1;SSC=16;GPV=1E0;SPV=2.4017E-2	GT:GQ:DP:RD:AD:FREQ:DP4	0/1:.:458:201:4:1.95%:196,5,0,4	0/1:.:1128:346:21:5.72%:342,4,0,21
    @@ -2766 +2766 @@
    -chr7	151970959	.	G	A	.	PASS	DP=3510;SS=1;SSC=62;GPV=2.4658E-14;SPV=5.4694E-7	GT:GQ:DP:RD:AD:FREQ:DP4	0/1:.:1216:1181:32:2.64%:420,761,8,24	0/1:.:2294:2281:13:0.57%:1799,482,4,9
    +chr7	151970959	.	G	A	.	PASS	DP=3510;SS=1;SSC=62;GPV=1E0;SPV=5.4694E-7	GT:GQ:DP:RD:AD:FREQ:DP4	0/1:.:1216:1181:32:2.64%:420,761,8,24	0/1:.:2294:2281:13:0.57%:1799,482,4,9

2. Add a variant named **freqQuot** in **comparePositions()** in *Somatic.java*
freqQuot= tumorFreq / normalFreq. 
In CASE 2E(Variant alleles match but difference significant), if freqQuot >= 10.0, we can treat this mutation as **Somatic**. For example tumorFreq=7.0% and normal=0.3%.

`diff -u after.vcf before.vcf`

    -chr2    29416481        .       T       C       .       PASS    DP=2030;SOMATIC;SS=2;SSC=198;GPV=1E0;SPV=1.3145E-20     GT:GQ:DP:RD:AD:FREQ:DP4 0/1:.:1078:1073:4:0.37%:859,214,3,1     0/1:.:952:877:74:7.78%:744,133,64,10
    +chr2    29416481        .       T       C       .       PASS    DP=2030;SS=1;SSC=198;GPV=1.5554E-24;SPV=1.3145E-20      GT:GQ:DP:RD:AD:FREQ:DP4 0/1:.:1078:1073:4:0.37%:859,214,3,1     0/1:.:952:877:74:7.78%:744,133,64,10
