#!/usr/bin/env python

import asyncio
from pathlib import Path
from playwright.sync_api import sync_playwright
import argparse

def get_args():
    parser = argparse.ArgumentParser(
        prog='pprint', # 程序名
        description='paged html 转 pdf'
    )
    parser.add_argument('-i', '--input', required=True,help="输入 html 文件",type=str)
    parser.add_argument('-o', '--output', required=True,help="输出 pdf 文件",type=str)
    args = parser.parse_args()
    return args


def convert(input,output):
    file = "file://" + str(Path(input).absolute())
    with sync_playwright() as p:
        browser_type = p.chromium
        browser = browser_type.launch()
        page = browser.new_page()
        page.goto(file)
        page.wait_for_selector(".pagedjs_pages")
        page.pdf(path=output,format="A4",outline = True,tagged = True)
        browser.close()
    return None

if __name__ == '__main__':
    args = get_args()
    convert(args.input,args.output)