Flaver Documentation

Description
Flaver: mining transcription factors in genome-wide transcriptome profiling data using weighted rank correlation statistics

Release
Source code	2.2.1	08/31/2022
Source code	2.2.0	06/14/2022
Binary	Windows, Mac OS, Linux

Install
Go the source folder, and run “make” command, a binary will produced in the “bin” folder

Requirement
Minimum 8GB memory, if you install from the source code, the g++ complier is also required.

Usage
flaver -s rank_set -i rank_data [-q q_score] -o output
Options
-m ranked gene sets
-i ranked genes with transcriptome data
-q fdr cut-off, optional
-o output, output file name

Example Data
rank_set file: human_rank_set.txt
rank_data file: skeletal-muscle_exp.txt
results: skeletal-muscle_output.txt 

Example
An example run should like: flaver -s human_rank_set.txt -i skeletal-muscle_exp.txt -o skeletal-muscle_output.txt
This command took two input files: human_rank_set.txt, skeletal-muscle_exp.txt. After finished run it will produce an output file named: skeletal-muscle_output.txt

References
Tinghua Huang, Xinmiao Huang, Binyu Wang, Hao He, Zhiqiang Du, Min Yao, and Xuejun Gao. Flaver: mining transcription factors in genome-wide transcriptome profiling data using weighted rank correlation statistics

Contact
Dr. Tinghua Huang, thua45@126.com
Dr. Min Yao, minyao@yangtzeu.edu.cn
