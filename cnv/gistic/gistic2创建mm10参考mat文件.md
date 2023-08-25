# 前言
- gistic2官方文件以及genepattern网站上面有的reference文件只有人类的hg16-hg38文件
- 需要创建不同物种的mat文件
- 采用模仿mat文件的策略：用matlab制作gistic mat参考文件

# 前期准备
- matlab
- gistic人类hg19.mat文件

## 人类hg19 mat参考文件解析
- 使用matlab打开hg19.mat文件
	- 工作区含有三个表格变量
		- cyto：用于记录cytoband的信息
			- chr：去除chr前缀的染色体名字（1:22，X，Y）；字符串
			- chrn：去除chr前缀的染色体名字，其中性染色体XY用23，24来代替；数字
			- name：cytoband名字（1p36.33等）；字符串
				- gisitc用name进行cytoband的注释
			- start：cytoband开始区域；数字
			- end：cytoband结束区域；数字
			- stain：吉姆萨染色法的结果，其中不同字符代表不同结果；字符串
		- rg：记录具体的转录本注释
			- refseq：转录本的refseq id；字符串
			- gene：对应基因具体名字；字符串
			- symb：对应基因symbol；字符串
			- locus_id：对应entrez id，基因的唯一标识
			- chr：含有chr前缀的染色体名字（chr1-chr22，chrX，chrY）；字符串
			- strand：基因所在的染色体的正链or负链（0,1）；数字
			- start：根据文件内相同基因的信息推断为该转录本的起始位置；数字
			- end：根据文件内相同基因的信息推断为该转录本的终止位置；数字
			- cds_start：编码dna序列起始位置；数字
			- cds_end：编码dna序列终止位置；数字
			- status：该注释的注释状态（例如Reviewed）；字符串
			- chrn：去除chr前缀的染色体名字，其中性染色体XY用23，24来代替；数字
		- rg_info：包含软件的metadata元数据，对于运行gistic2来说并不需要，但是需要用来编译更新以后的gp_gistic2_from_seg可执行脚本
			- Field一列，Value一列
				- assembly：GRCh37/hg19
				- source：TCGA
				- path：/xchip/gistic/variables/20120227_UCSC_dump_hg19
				- url：fftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/
				- script：'/xchip/gistic/variables/20120227_UCSC_dump/hg19/make_rg_20120227.m'
				- makeversion：“”
				- maker：schum
				- date：27-Feb-2012

# 小鼠mm10 mat文件制作与重新编译
## 小鼠mm10 mat参考文件制作（GRCm38）
- cyto表格的制作
	- 从UCSC网站中下载cytoband table（UCSC table browser或者UCSC goldenPath archive）：[http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/cytoBand.txt.gz](http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/cytoBand.txt.gz)
- rg表格制作
	- 方法一：下载基因组版本的具体基因注释文件：refGene、Ensembl、GENCODE中的gtf文件（未使用过）
	- 方法二
		- [Table Browser (ucsc.edu) | 使用UCSC table browser做](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1682280320_PyXA25KMEfQIMXupGWI6kZ0Os7bu&clade=mammal&org=Mouse&db=mm10&hgta_group=genes&hgta_track=refSeqComposite&hgta_table=refGene&hgta_regionType=genome&position=chr12%3A56%2C694%2C976-56%2C714%2C605&hgta_outputType=primaryTable&hgta_outFileName=result2.txt)
			- 参数如下
				- clase：mammal
				- genome：mouse
				- assembly：GRCm38/mm10
				- group：genes and gene predictions
				- track：NCBI refseq
				- table：refseq all
			- 结果：获得小鼠中mm10所有的refseq的mrna以及ncrna的信息及其对应的cds信息、基因信息（主要信息）
		- 使用biomarts进行转id，获取剩余信息
- rg_info表格制作
- mat文件制作

```R
# this script is for creating csv files for RG and RG_INFO
# necessary packages
packages <- c("dplyr","stringr","tidyverse","biomaRt","rtracklayer")
for (i in packages) {
  if (!require(i, character.only = TRUE)) {
    stop(paste("package: ", i, " is not exits"))
  }
}

# biomart for converting the refid
listMarts(host="https://nov2020.archive.ensembl.org/")
mm10_mart <- useMart("ensembl",host="https://nov2020.archive.ensembl.org/")
listDatasets(mm10_mart)
mm10_mart <- useDataset("mmusculus_gene_ensembl", mm10_mart)
listFilters(mm10_mart)
View(listAttributes(mm10_mart))

refseq <- read.table("C:/Users/Administrator/Desktop/result (1).txt",header=TRUE) %>%
  dplyr::select(name,cdsStart,cdsEnd) %>%
  dplyr::filter(stringr::str_detect(name,pattern="^NM")) %>%
  dplyr::mutate(name=stringr::str_remove_all(name,pattern="\\..*")) %>%
  dplyr::distinct() %>%
  dplyr::left_join(biomaRt::getBM(attributes=c("refseq_mrna","description","external_gene_name","chromosome_name",
                                               "transcript_start","transcript_end","entrezgene_id",
                                               "strand"),filter="refseq_mrna",
                                  values=.$name[grep("^NM",.$name)] %>% unlist() %>% as.character(),mart=mm10_mart),
                   by=c("name"="refseq_mrna")) %>%
  dplyr::bind_rows(read.table("C:/Users/Administrator/Desktop/result (1).txt",header=TRUE) %>%
                     dplyr::select(name,cdsStart,cdsEnd) %>%
                     dplyr::filter(stringr::str_detect(name,pattern="^NR")) %>%
                     dplyr::mutate(name=stringr::str_remove_all(name,pattern="\\..*")) %>%
                     dplyr::distinct() %>%
                     dplyr::left_join(biomaRt::getBM(attributes=c("refseq_ncrna","description","external_gene_name","chromosome_name",
                                                                  "transcript_start","transcript_end","entrezgene_id",
                                                                  "strand"),filter="refseq_ncrna",
                                                     values=.$name[grep("^NR",.$name)] %>% unlist() %>% as.character(),mart=mm10_mart),
                                      by=c("name"="refseq_ncrna"))) %>%
  dplyr::filter(!is.na(description)) %>%
  dplyr::filter(!is.na(external_gene_name)) %>%
  dplyr::filter(chromosome_name %in% c(1:19,"X","Y")) %>%
  dplyr::filter(!is.na(entrezgene_id)) %>%
  dplyr::rename(refseq=name,
                gene=description,
                symb=external_gene_name,
                locus_id=entrezgene_id,
                chr=chromosome_name,
                strand=strand,
                start=transcript_start,
                end=transcript_end,
                cds_start=cdsStart,
                cds_end=cdsEnd) %>%
  dplyr::mutate(chrn=as.numeric(unlist(lapply(chr,function(x){ifelse(x=="X",20,ifelse(x=="Y",21,x))}))),
                chr=paste("chr",chr,sep = ""),
                status="Reviewed",
                strand=unlist(lapply(strand,function(x){ifelse(x==-1,0,1)}))) %>%
  dplyr::select(refseq,gene,symb,locus_id,chr,strand,start,end,cds_start,cds_end,status,chrn)
                            
write.table(refseq,"C:/Users/Administrator/Desktop/rg.txt",sep="\t",quote = FALSE,row.names = FALSE)    

```

```matlab
### first read cyto table
# ps: because chr column is a mixture of numeric and character, so we need to set options
opt=detectImportOptions("C:\Users\Administrator\Desktop\cytoBand.csv")
opt=setvartype(opt,{'chr'},'char')
opt=setvartype(opt,{'start'},'int32')
opt=setvartype(opt,{'xEnd'},'int32')
cyto=readtable("C:\Users\Administrator\Desktop\cytoBand.csv",opt)
cyto = sortrows(cyto,"chr","ascend")
cyto = table2struct(cyto)

### then read eg table
opt=detectImportOptions("C:\Users\Administrator\Desktop\rg.txt")
opt=setvartype(opt,{'chr'},'char')
opt=setvartype(opt,{'start'},'int32')
opt=setvartype(opt,{'xEnd'},'int32')
opt=setvartype(opt,{'cds_start'},'int32')
opt=setvartype(opt,{'cds_end'},'int32')
rg=readtable("C:\Users\Administrator\Desktop\rg.txt",opt) 
rg = table2struct(rg)

### then read rg_info table
rg_info=struct(assembly=char("GRCm38/mm10"),source=char("UCSC"),path=char("https://genome.ucsc.edu/"),url=char(""),script=char(""),make_version=char(""),maker=char("champeil"),date=char("21-Aug-2023"))

# save mat file
save C:\Users\Administrator\Desktop\mm10.mat cyto rg rg_info
```

```python
# 由于在matlab中struct结构的名字不能改成end（跟读取函数的参数有冲突），所以需要从python中读取并修改
import scipy.io
import numpy as np
mat_data=scipy.io.loadmat("/home/laojp/software/gistic_mm10/refgenefiles/mm10.mat")
mat_data['cyto'].dtype.names=('chr', 'chrn', 'name', 'start', 'end', 'stain')
mat_data['rg'].dtype.names=('refseq','gene','symb','locus_id','chr','strand','start','end','cds_start','cds_end','status','chrn')
scipy.io.savemat("/home/laojp/software/gistic_mm10/refgenefiles/mm10.mat", mat_data)

```
## gp_gistic2_from_seg脚本重新编译
- 前言
	- 如果使用非人类的gistic参考mat文件，则需要重新编译gp_gistic2_from_seg脚本，因为某些gistic2源代码文件（source/RefGeneInfo.m）含有人类基因组的硬编码染色体
- 方法：使用matlab重新编译gp_gistic2_from_seg
- 过程
	- 修改source/RefGeneInfo.m文件
	- 重新编译gp_gistic2_from_seg 模块
```shell
1. 打开source/RefGeneInfo.m
	1. 定位到15行，将nchr改成目标物种含有的染色体数目（小鼠含有21条，所以变成21）
		nchr = 21;
	2. 定位到29行，添加、修改对应的染色体名字（小鼠变成1:19，X，Y）
		RGI.chr.symb = {'1','2','3','4','5','6','7','8','9','10',...'11','12','13','14','15','16','17','18',...'19','X','Y'};
	3. 定位到40-42行，以匹配染色体数量
		RGI.txt2num('20') = 20;   
		RGI.txt2num('21') = 21;    
		RGI.chr.autosomal = (1:21)<20; 
```
- 关于重新编译的方法
	- linux中安装matlab（gistic2使用的是2014a）
	- 调用matlab环境
	- 使用matlab中mcc编译source/gp_gistic2_from_seg.m脚本，得到新的gp_gistic2_from_seg可执行文件
	- 将新的可执行文件替换掉gistic2脚本中的文件
