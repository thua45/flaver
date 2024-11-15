Flaver Documentation

Description
Testing the significance of ranked gene-sets in genome-wide transcriptome profiling data using weighted rank correlation statistics

Release
Source code	1.0.31	08/18/2024
Source code	1.0.2024	06/18/2023
Source code	1.0.1	08/31/2022
Source code	1.0.0	06/14/2022
Binary	Windows, Mac OS, Linux

Install
Go the source folder and run “g++ main.cpp -std=c++11 -lpthread -o flaver” command, a binary named “flaver” will produced in the folder.
Under windows OS try “g++ main.cpp -std=c++11 -static -lpthread -o flaver.exe” instead.

Requirement
Minimum 2GB memory, if you install from the source code, the g++ complier is also required.

Usage
flaver -s gene_set -i gene_list [-g tscore] [-r adjust] [-w weight] [-p auto_p] [-t thread] [-o output] [-d outdata]
Options
-s ranked gene-sets.
-i ranked gene-list.
-g gene-set score, less than this value will be omitted in the analysis, default 10.
-r two letter parameter for data adjust (first letter for gene-set, second letter for gene-list), 1 for as is, 2 for absolute value, and 3 for reverse value (x * -1.0), default "11". Example of "12" will be take the gene-set value as is and the absolute values for gene-list.
-w weighting method, 1 by geometric mean of gene-set and gene-list, 2 by gene-list, 3 by gene-set, 4 by geometric mean of gene-set and gene-list (ori), 5 by gene-list (ori), 6 by gene-set (ori), 7 by mixed-density, 8 by gene-list-density, 9 by gene-set-density, default 1.
-p auto p, default is off.
-t multithreading, default 1. Usually set to the number of CPU cores.
-o output, output file name.
-d output data file name.

Example Data
ranked gene-set file: Grit-set
ranked gene-list file: RNA55T-lists
results: Grit-set-2-RNA55T-lists

Example
An example run should like: ./flaver -s grit-human-v100-j20-2.bed -i duodenum_gapsig-e6.txt -g 1.3 -w 3 -t 8 -o duodenum_ output.txt
This command took two input files: grit-human-v100-j20-2.bed, duodenum_gapsig-e6.txt. After finished run it will produce an output file nameded duodenum_ output.txt.

Benchmark
The benchmark result for human transcriptome with Flaver is available for download.

References
Min Yao, Hao He, Binyu Wang, Xinmiao Huang, Sunli Zheng, Jianwu Wang, Xuejun Gao, Tinghua Huang. Testing the significance of ranked gene-sets in genome-wide transcriptome profiling data using weighted rank correlation statistics.

Contact
Dr. Tinghua Huang, thua45@126.com
Dr. Min Yao, minyao@yangtzeu.edu.cn
