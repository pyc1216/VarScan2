# VarScan2
version=v2.3.9.b

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
    

2. Fix a bug: When **normalReads2** = 1 < **minReads2** = 2, this may make **diffPvalue** > *0.05*. Set **normalReads2** = 0, if normalAllele2.equals(refBase).

#### my version:
    chr7	55242464	.	AGGAATTAAGAGAAGC	A	.	PASS	DP=1819; SOMATIC;SS=2;SSC=6;GPV=1E0;SPV=2.3063E-1	GT:GQ:DP:RD:AD:FREQ:DP4	0/0:.:305:304:0:0%:172,132,1,0	0/1:.:1514:1506:8:0.53%:692,814,6,2
    
#### origin version:
    chr7    55242464	.	AGGAATTAAGAGAAGC	A	.	PASS	DP=1819;SS=0;SSC=0;GPV=1E0;SPV=5.3897E-1	GT:GQ:DP:RD:AD:FREQ:DP4	0/0:.:305:304:1:0.33%:172,132,1,0	0/1:.:1514:1506:8:0.53%:692,814,6,2
    
