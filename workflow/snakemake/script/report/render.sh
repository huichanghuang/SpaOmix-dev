#!/usr/bin/env bash

# 一次性生成 html、pdf、word报告
quarto render

# word 报告添加首页
python script/add_docx.py -i report.docx -o report.docx

# html_paged 转 pdf
