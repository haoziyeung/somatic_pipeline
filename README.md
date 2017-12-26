# somatic_pipeline
肿瘤体细胞突变检测流程（组织、重复DNA均可，需要有对照）
<pre>
----------------------------------------------------------
        {体细胞突变检测流程}
特点：
1)使用sentieon TNscope流程，准确、高效、快速
2)关联cosmic、www.cancergenomeinterpreter.org靶药数据库
3)分析SV，常见的融合基因一目了然
4)分析CNV，可检测拷贝数变异


输入：fastq

输出：somatic_mut.xlsx/tsv,SV,CNV,QC metrics...

使用方法：
        perl somatic.sentieon.fq.pipeline.pl -t <tumor name> -t1 <tumor R1> -t2 <tumor R2> -n <normal name> -n1 <normal R1> -n2 <normal R2> -b <interval.bed>

默认提交到cn-medium分区，如果没有资源了，你可以指定作业提交到cn-fast，如下所示：
        perl somatic.sentieon.fq.pipeline.pl -t <tumor name> -t1 <tumor R1> -t2 <tumor R2> -n <normal name> -n1 <normal R1> -n2 <normal R2> -b <interval.bed> -part cn-fast

Last Updated:2017/12/25

如果使用过程中发生任何问题，请联系yanghao@eulertechnology.com

---------------------------------------------------------
</pre>
