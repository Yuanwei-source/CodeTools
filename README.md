<<<<<<< HEAD
# Annotation_protein_filter.py 使用说明

## 简介

本脚本用于对蛋白质序列（FASTA格式）进行多步过滤和处理，包括：
- 移除含有非标准氨基酸的序列
- 按长度过滤序列
- 使用 CD-HIT 进行聚类去冗余
- 重命名重复的序列ID
- 生成处理报告

## 依赖

- Python 3
- [Biopython](https://biopython.org/)
- [seqkit](https://bioinf.shenwei.me/seqkit/)
- [CD-HIT](http://weizhongli-lab.org/cd-hit/)

## 参数说明

| 参数 | 说明 | 默认值 |
|------|------|--------|
| `-i, --input` | 输入FASTA文件 | input.faa |
| `-o, --output` | 输出FASTA文件 | processed.faa |
| `-r, --report` | 记录被移除序列ID的文件 | removed.txt |
| `-t, --threads` | 线程数 | CPU核心数 |
| `-c, --chunk` | 每个处理块的序列数 | 50000 |
| `--keep-uo` | 保留U和O氨基酸 | 无 |
| `--min-len` | 最小序列长度 | 100 |
| `--max-len` | 最大序列长度 | 1500 |
| `--identity` | CD-HIT聚类相似性阈值 | 0.90 |
| `--skip-length` | 跳过长度过滤 | 无 |
| `--skip-cdhit` | 跳过CD-HIT聚类 | 无 |

## 使用示例

```bash
python Annotation_protein_filter.py -i input.faa -o output.faa -r removed.txt --min-len 120 --max-len 2000 --identity 0.95 --threads 8
```

## 输出文件

- 过滤后的FASTA文件（`-o`指定）
- 被移除序列的ID及原因（`-r`指定）
- CD-HIT聚类结果（`output.faa.clstr`，如未跳过CD-HIT）

## 注意事项

- 需提前安装 seqkit 和 cd-hit，并确保可在命令行调用。
- 输入文件需为标准FASTA格式。
- 支持多线程加速处理。

=======
# CodeTools
some toosl for bioinformatics
>>>>>>> 519d5ff76628cc2d6be0f6bedd87d966c2002d1d
