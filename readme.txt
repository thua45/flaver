Flaver-2.0 Documentation

Description
Testing the significance of ranked gene-sets in genome-wide transcriptome profiling data using weighted rank correlation statistics

Homepage
http://www.thua45.cn/grit-2.0

Release
Binary-2.0	Windows	01/12/2025	http://www.thua45.cn/grit-2.0/programs/flaver-2.0/flaver-2.0.zip
Source code	1.0.31	08/18/2024
Source code	1.0.2024	06/18/2023
Source code	1.0.1	08/31/2022
Source code	1.0.0	06/14/2022
Binary-1.0	Windows, Mac OS, Linux

Install
Go the source folder and run "g++ main.cpp -std=c++11 -lpthread -o flaver" command, a binary named "flaver" will produced in the folder.
Under windows OS try "g++ main.cpp -std=c++11 -static -lpthread -o flaver.exe" instead.

Requirement
Minimum 2GB memory, if you install from the source code, the g++ complier is also required.

Usage
flaver -s gene_set -i gene_list [-t thread] [-o output]

Options
-s ranked gene-sets.
-i ranked gene-list.
-t multithreading, default 1. Usually set to the number of CPU cores.
-o output, output file name.

Example Data
ranked gene-set: Pig10.2-FLAVER.bed	http://www.thua45.cn/grit-2.0/programs/flaver-2.0/Pig10.2-FLAVER.bed
ranked gene-list: gene_list.txt	http://www.thua45.cn/grit-2.0/programs/flaver-2.0/gene_list.txt

Example
An example run should like: flaver-2.0 -s Pig10.2-FLAVER.bed -i gene_list.txt -t 8 -o output.txt.
This command takes two input files: Pig10.2-FLAVER.bed and gene_list.txt. It will produce an output file named output.txt. Full instruction is avaiavle in Flaver's user manual.

References
Min Yao, Hao He, Binyu Wang, Xinmiao Huang, Sunli Zheng, Jianwu Wang, Xuejun Gao, Tinghua Huang. Testing the significance of ranked gene-sets in genome-wide transcriptome profiling data using weighted rank correlation statistics, Current Genomics, 2024.
Huang Tinghua, Niu Siqi, Zhang Fanghong, Wang Binyu, Wang Jianwu, Liu Guoping, Yao Min, Correlating gene expression levels with transcription factor binding sites facilitates identification of key transcription factors from transcriptome data, Frontiers in Genetics, 2024.

Contact
Dr. Tinghua Huang, thua45@126.com
Dr. Min Yao, minyao@yangtzeu.edu.cn
