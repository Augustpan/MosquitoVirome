准备四个文件：
    1. query.fa 自己的病毒蛋白序列（RdRp、DNA病毒复制酶）
    2. ref.fa 画树的参考序列
    3. manual_assignment.csv 参见示例文件
    4. virus_manifest.txt 所有自己的病毒的名字

先跑：sh blast_search.sh query.fa ref.f，然后跑R

输出两个文件：
    1. OUTPUT_HOST_TABLE.csv 自己的病毒的宿主注释
    2. OUTPUT_HOST_TABLE_WITH_REF.csv 自己的病毒+参考序列的注释（用来画树）
