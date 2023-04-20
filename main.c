#include <stdio.h>
#include <unistd.h>
#include<string.h>

#ifndef PACKAGE_VERSION
#define PACKAGE_VERSION "0.0.1"
#endif

int kmer_count(int argc, char *argv[]);
int kmer_match(int argc, char *argv[]);
int kmer_dump(int argc,char *argv[]);

int kmer_set_group(int argc,char *argv[]);
int kmer_set_intersect(int argc,char *argv[]);
int kmer_set_union(int argc,char *argv[]);
int kmer_set_diff(int argc,char *argv[]);

static int usage()
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Program: tod (Tools for kmer analysis)\n");
	fprintf(stderr, "Version: %s\n\n", PACKAGE_VERSION);
	fprintf(stderr, "Usage:   tod <command> [options]\n\n");
	fprintf(stderr, "Command: count            count kmers \n");
	fprintf(stderr, "         dump             dump kmers in a format :k-mer counts [counts...]\n");
	fprintf(stderr, "         group            group multiple kmer files\n");
	fprintf(stderr, "         intersect        intersection of multiple kmer files\n");
	fprintf(stderr, "         union            union of multiple kmer files\n");
	fprintf(stderr, "         diff             difference of multiple kmer files,use the first file as a control\n");
	fprintf(stderr, "         map              map sequences to the target genome\n");
	fprintf(stderr, "\n");
	return 1;
}

/*
tod count -k 31 -t 4 [-h out.hist] file.fa > file.kmc 
    切kmer 输出的kmc文件用dump命令查看
    支持多线程
    最大支持31mer 计数最大支持16383
    
tod paste -f 'cnt[1]>1 and cnt[2]=0  and cnt[1:3]>10 and number[1,2,3]>=2' file1.kmc file2.kmc file3.kmc 
    kmer合并 可完成样品的和差并等运算
    cnt为一组样品的kmer计数
    number为满足kmer频次要求的样品数,值为小数时表示出现比例
    支持[1:5,9]这种索引,闭区间 
tod paste -l A,A,B -f 'cnt[A]>=2 and cnt[B]<20' file1.kmc file2.kmc file3.kmc 
    指定样品名称，样品名称相同的在输出文件中合并为一列
    
tod map -t 4 -d target.fa  <qry.fa|qry.kmc|kmer.list>
    查找指定kmer/fa在target.fa中的位置
    当输入为kmer文件时，输出结果kmer tseqname1 count position信息文件
    	ACG tseq1 3 7,23,100
    	ACG tseq2 4 5,67,89,105,233
    当输入为fa文件时，输出结果为qseqNameList tseqname pos信息文件
        qseq1,qseq2,qseq3 tseq1 100 2100

tod intersect file1 file2 [-c 1,2] [-C 100,200] -o collapse|add|sub|min|max|1/2 
	输出file1 file2共有kmer 
	计数结果由输出选项控制 默认为空表示不输出
	其他输出选项collapse表示折叠输出全部计数 add/sub/min/max表示输出两个kmer计数的和/差/较小/较大值 1/2表示输出1/2中的kmer计数
	
tod union file1 file2 [-c 1,2] [-C 100,200] -o collapse|add|sub|min|max|1/2 
    输出file1 file2 至少一个文件中存在的kmer
	
tod diff file1 file2  [-c 1,2] [-C 100,200] -o collapse|add|sub|min|max|1/2 
    输出file1中有file2中没有的kmer
	
tod group file1 file2 file3 file4 [-c 1,2,2,1] [-C 100,200,200,100] [-n 2]  [-N 3] -o collapse|add|sub|min|max|[id] 
	输出以file1/2/3/4作为一个组，限定该组中个体数为2-3个的kmer


联合使用：
tod intersect <(tod group file1 file2 file3 -n2) <(tod group file4 file5 file6 -N1) 即在file1/2/3中至少两个文件存在，在file4/5/6中至多1个文件存在
tod dump file.kmc
	以tab格式 查看计数文件
    
*/
int main(int argc, char *argv[])
{
	if (argc < 2) return usage();
	if (strcmp(argv[1], "count") == 0) return kmer_count(argc-1, argv+1);
	else if (strcmp(argv[1], "dump") == 0) return kmer_dump(argc-1, argv+1);
	else if (strcmp(argv[1], "group") == 0) return kmer_set_group(argc-1, argv+1);
        else if (strcmp(argv[1], "intersect") == 0) return kmer_set_intersect(argc-1, argv+1);
        else if (strcmp(argv[1], "union") == 0) return kmer_set_union(argc-1, argv+1);
        else if (strcmp(argv[1], "diff") == 0) return kmer_set_diff(argc-1, argv+1);
	else if (strcmp(argv[1], "map") == 0) return kmer_match(argc-1, argv+1);
	else {
		fprintf(stderr, "[main] unrecognized command '%s'\n", argv[1]);
		return 1;
	}
	return 0;	
}

