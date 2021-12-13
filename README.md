# ANNOVAR tools

该包是配合[ANNOVAR](https://annovar.openbioinformatics.org/en/latest/) 对相关文件进行处理的工具集

## vcf

> 将ANNOVAR输入文件,AVinput文件转换为VCF4.0格式

- -i/--avinput: 输入文件,AVinput格式
- -r/--reference: 参考基因组文件如hg19.fa,请事先使用samtools faidex构建索引
- -o/--vcf: 输出的VCF文件

```shell
zhuy@ubuntu:/projects/example$ annovar_tools.py vcf -i test.avinput -r hg19.fa -o test.vcf
```

## split

> 将ANNOVAR注释结果按照Gene拆分为多行，一次只能拆分一种gene-based数据库

- -i/--avoutput: 输入文件,ANNVOAR注释结果，如test.hg19_multianno.txt
- -r/--refgenes: 运行ANNOVAR所用的注释refGene文件，如hg19_ensGene.txt,hg19_knownGene.txt,hg19_refGeneWithVer.txt
- -g/gene_db: 所需拆分的gene-based数据库, 如refGeneWithVer,即ANNVOAR g参数对应的数据库名称,在ANNOVAR结果中体现为"Func.refGeneWithVer, Gene.refGeneWithVer, GeneDetail.refGeneWithVer, ExonicFunc.refGeneWithVer,
  AAChange.refGeneWithVer"
- -o/--outfile: 输出的拆分后的文件

```shell
zhuy@ubuntu:/projects/example$ annovar_tools.py split \
                               -i test.hg19_multianno.txt \
                               -r hg19_refGeneWithVer.txt \
                               -r hg19_ensGene.txt \
                               -r hg19_knownGene.txt \
                               -g refGeneWithVer
                               -o test.hg19_multianno.refGeneWithVer.txt
```