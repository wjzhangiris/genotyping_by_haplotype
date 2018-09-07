# genotyping_by_haplotype
using methylation haplotype to infer CNV genotype

## 利用下面脚本推测目标区域CNV的基因型
genotype_by_target.py


## 程序用法
python genotype_by_target.py target.config


## target.config文件
该文件是json格式，例子文件参照config/target.config
{
"raw_hap_plus":"",
"raw_hap_minus":"",
"villus_input":"",
"pregnant_plasma_input":"",
"villus_sample":"",
"pregnant_plasma_sample":"",
"output":"",
"target_hapid":""
}
