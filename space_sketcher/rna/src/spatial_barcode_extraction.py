import os
import sys
import argparse
import dnaio
import pandas as pd
import numpy as np
from itertools import takewhile
from typing import Tuple, Dict, Set, List
import polars as pl
from space_sketcher.tools.utils import add_log


CHUNK_SIZE = 5_000_000    # 每个chunk处理5M reads

def read_whitelist_files(file_path: str) -> Set[str]:
    """优化版白名单读取，使用生成器表达式减少内存使用"""
    with open(file_path, 'r') as f:
        return {line.split()[0] for line in f if line.strip()}

def get_whitelist(cb_file: str, chem: str, params: str = None) -> Set[str]:
    """使用numpy加速白名单处理"""
    base_wl = read_whitelist_files(cb_file)
    if chem == "10X":
        return base_wl
    elif chem == "leader_v1":
        arr = np.array(list(base_wl))
        return {'_'.join(pair) for pair in np.array(np.meshgrid(arr, arr)).T.reshape(-1, 2)}
    elif chem == "other":
        if not params:
            raise ValueError("Missing mapparams while rnachemistry is other!")
        
        parts = params.split()
        if "--soloCBwhitelist" not in parts:
            raise ValueError("Missing --soloCBwhitelist in mapparams!")
        
        wlindex = parts.index("--soloCBwhitelist")
        cb_files = list(takewhile(lambda x: not x.startswith("--"), parts[wlindex + 1:]))
        
        if len(cb_files) > 2:
            raise ValueError("Does not support more than 2 cell barcode whitelists!")

        if "CB_UMI_Simple" in params:
            return read_whitelist_files(cb_files[0])
        elif "CB_UMI_Complex" in params:
            barcode_lists = [read_whitelist_files(f) for f in cb_files]
            arrays = [np.array(list(b)) for b in barcode_lists]
            grids = np.meshgrid(*arrays, indexing='ij')
            return {'_'.join(combo) for combo in np.stack(grids, -1).reshape(-1, len(arrays))}
        else:
            raise ValueError("Invalid mapparams format")
    else:
        raise ValueError(f"Invalid rnachemistry: {chem}")

def check_rnachemistry(chem_type: str, mapparams: str = None) -> Dict:
    """验证RNA化学类型并提取位置信息"""
    pos_info = {
        "cbstart": [],
        "cblen": [],
        "umistart": 0,
        "umilen": 0
    }
    
    try:
        if chem_type == "10X":
            pos_info.update({
                "cbstart": [0],
                "cblen": [16],
                "umistart": 16,
                "umilen": 12
            })
        elif chem_type == "leader_v1":
            pos_info.update({
                "cbstart": [0, 10],
                "cblen": [10, 10],
                "umistart": 20,
                "umilen": 10
            })
        elif chem_type == "other":
            if not mapparams:
                raise ValueError("Missing mapparams while rnachemistry is other!")
            
            parts = mapparams.split()
            if not validate_mapparams(parts):
                raise ValueError("Invalid mapparams structure")
                
            if "CB_UMI_Simple" in mapparams:
                process_simple_case(parts, pos_info)
            elif "CB_UMI_Complex" in mapparams:
                process_complex_case(parts, pos_info)
            else:
                raise ValueError("Unsupported barcode type in mapparams")
        else:
            raise ValueError(f"Unsupported library type: {chem_type}")
    except Exception as e:
        logger.error(f"Error in check_rnachemistry: {str(e)}")
        sys.exit(1)
        
    return pos_info

def validate_mapparams(parts: List[str]) -> bool:
    """验证mapparams结构"""
    required_params = {
        "CB_UMI_Simple": ["--soloCBstart", "--soloCBlen", "--soloUMIstart", "--soloUMIlen"],
        "CB_UMI_Complex": ["--soloCBposition", "--soloUMIposition"]
    }
    
    for chem, params in required_params.items():
        if chem in ' '.join(parts):
            return all(p in parts for p in params)
    return False

def process_simple_case(parts: List[str], pos_info: Dict) -> None:
    """处理简单CB_UMI情况"""
    pos_info["cbstart"] = [int(parts[parts.index("--soloCBstart") + 1]) - 1]
    pos_info["cblen"] = [int(parts[parts.index("--soloCBlen") + 1])]
    pos_info["umistart"] = int(parts[parts.index("--soloUMIstart") + 1]) - 1
    pos_info["umilen"] = int(parts[parts.index("--soloUMIlen") + 1])

def process_complex_case(parts: List[str], pos_info: Dict) -> None:
    """处理复杂CB_UMI情况"""
    cb_index = parts.index("--soloCBposition")
    cb_values = list(takewhile(lambda x: not x.startswith("--"), parts[cb_index + 1:]))
    
    if len(cb_values) > 2:
        raise ValueError("Does not support more than 2 part of cell barcode")
    
    for cb_val in cb_values:
        _, start, _, end = cb_val.split('_')
        pos_info["cbstart"].append(int(start))
        pos_info["cblen"].append(int(end) - int(start) + 1)
        
    pos_info["umistart"] = int(parts[parts.index("--soloUMIstart") + 1]) - 1
    pos_info["umilen"] = int(parts[parts.index("--soloUMIlen") + 1])


def hamming_match_position_levenshtein(seq: str, pattern: str, max_mismatches: int = 1) -> int:
    L = len(pattern)
    N = len(seq)
    if N < L:
        return -1
    for i in range(N - L + 1):
        mismatches = 0
        for j in range(L):
            if seq[i+j] != pattern[j]:
                mismatches += 1
                if mismatches > max_mismatches:
                    break
        if mismatches <= max_mismatches:
            return i
    return -1

def find_linker(seq: str, pattern: str, max_mismatches: int = 1) -> Tuple[int, bool]:
    # 精确匹配
    pos = seq.find(pattern)
    if pos != -1:
        return pos

    # 汉明距离匹配
    pos = hamming_match_position_levenshtein(seq, pattern, max_mismatches)
    if pos != -1:
        return pos

    return -1

def auto_detect_sb(seq: str, linkers: list) -> Tuple[str, bool]:
    """
    Return (specific binding sequence, is_perfect_match)
    is_perfect_match = True only if all involved linkers matched exactly.
    """

    if len(linkers) == 2:
        linker1_pos = find_linker(seq, linkers[0])
        if linker1_pos == -1:
            return ""
        search_start = linker1_pos + len(linkers[0])
        linker2_pos = find_linker(seq[search_start:], linkers[1])
        if linker2_pos == -1:
            return ""
        linker2_pos += search_start

        try:
            sb1 = seq[linker1_pos - 6:linker1_pos]
            sb2 = seq[linker1_pos + len(linkers[0]):linker1_pos + len(linkers[0]) + 8]
            sb3 = seq[linker2_pos + len(linkers[1]):linker2_pos + len(linkers[1]) + 8]
            return sb1 + sb2 + sb3
        except IndexError:
            return ""

    elif len(linkers) == 1:
        linker_pos, perfect = find_linker(seq, linkers[0])
        if linker_pos != -1 and linker_pos >= 30:
            return seq[linker_pos - 30:linker_pos], perfect
        else:
            return ""

    return ""

@add_log
def process_fastq_singlethread(
    r1_path: str,
    r2_path: str,
    oligochip: str,
    cb_positions: List[Tuple[int, int]],
    umi_pos: Tuple[int, int],
    cb_whitelist: Set[str],
    sb_whitelist: Set[str] = None
):
    """
    单进程顺序处理 FASTQ 文件，返回符合条件的 (CB, UMI, SB, count=1) 元组。
    """
    results = []
    umi_start, umi_len = umi_pos
    check_sb = oligochip == "LD"
    linkerinfo = {
        "LD": ['TCTTCAGCGTTCCCGAGATCGGACGATCATGGG', 'CAAGTATGCAGCGCGCTCAAGCACGTGGAT'],
        "GM": ['TCTTGTGACTACAGCACCCTCGACTCTCGC']
    }
    if oligochip not in linkerinfo:
        raise ValueError("Invalid oligochip")

    linkers = linkerinfo[oligochip]
    total = cbmatched = 0
    with dnaio.open(r1_path, r2_path, mode="r") as reader:
        for r1, r2 in reader:
            total += 1
            if total % 1_000_000 == 0:
                print(f"Processed {total:,} read pairs...")
            cb_parts = [r1.sequence[start:start+length] for start, length in cb_positions]
            cb = "_".join(cb_parts) if len(cb_parts) > 1 else cb_parts[0]

            if cb not in cb_whitelist:
                continue
            cbmatched += 1

            umi = r1.sequence[umi_start:umi_start+umi_len]
            sb = auto_detect_sb(r2.sequence, linkers)

            if check_sb and sb_whitelist and sb not in sb_whitelist:
                continue

            results.append(f"{cb},{umi},{sb},1")
          
    csv_text = "Cell_Barcode,UMI,Spatial_Barcode,Read_Count\n" + "\n".join(results)

    # 使用 Polars 从字符串读取 CSV
    results = pl.read_csv(
        source=csv_text.encode(),
        separator=","
    )
    results = results.group_by(["Cell_Barcode", "UMI", "Spatial_Barcode"]).agg(
        pl.col("Read_Count").sum()
    )
    
    print(f"Total reads processed: {total:,}, Cell barcode matched reads: {cbmatched:,}")
    return results,total, cbmatched

@add_log
def process_fastq_multithread(
    r1_path: str,
    r2_path: str,
    oligochip: str,
    cb_positions: List[Tuple[int, int]],
    umi_pos: Tuple[int, int],
    cb_whitelist: Set[str],
    sb_whitelist: Set[str] = None,
    cores: int = 4,
    outdir: str = None
):
    from cutadapt.runners import make_runner
    from cutadapt.files import InputPaths
    from cutadapt.pipeline import Pipeline
    from cutadapt.utils import Progress
    input_paths = InputPaths(
        *[
            r1_path,
            r2_path,
        ],
    )
    class NullOutputFiles:
        def proxy_files(self):
            return []

        def binary_files(self):
            return []

    check_sb = oligochip == "LD"
    linkerinfo = {
        "LD": ['TCTTCAGCGTTCCCGAGATCGGACGATCATGGG', 'CAAGTATGCAGCGCGCTCAAGCACGTGGAT'],
        "GM": ['TCTTGTGACTACAGCACCCTCGACTCTCGC']
    }
    if oligochip not in linkerinfo:
        raise ValueError("Invalid oligochip")

    linkers = linkerinfo[oligochip]
    if outdir is None:
        raise ValueError("Output directory must be specified for multithreaded processing")
    temp_dir = os.path.join(outdir, "_multithread_temp")

    # 如果目录已存在，先删除
    if os.path.exists(temp_dir):
        import shutil
        shutil.rmtree(temp_dir)
    os.makedirs(temp_dir)

    class SimplePipeline(Pipeline):
        def __init__(self):
            self.paired = True
            self._modifiers = []
            self._steps = []
            self.cb_positions = cb_positions
            self.umi_pos = umi_pos
            self.cb_whitelist = cb_whitelist
            self.sb_whitelist = sb_whitelist
            self.linkers = linkers
            self.check_sb = check_sb
            self.outdir = outdir

        def process_reads(self, infiles, progress=None):
            pid = os.getpid()
            self._infiles = infiles
            self._reader = infiles.open()
            results = []
            n = perfect_match_count = cbmatched = 0
            for reads in self._reader:
                n += 1
                if n % 10000 == 0 and progress is not None:
                    progress.update(10000)
                r1, r2 = reads
                cb_parts = [r1.sequence[start:start+length] for start, length in self.cb_positions]
                cb = "_".join(cb_parts) if len(cb_parts) > 1 else cb_parts[0]

                if cb not in cb_whitelist:
                    continue
                cbmatched += 1

                umi_start, umi_len = self.umi_pos
                umi = r1.sequence[umi_start:umi_start+umi_len]
                sb = auto_detect_sb(r2.sequence, linkers)
                if check_sb and sb_whitelist and sb not in sb_whitelist:
                    continue

                results.append(f"{cb},{umi},{sb},1\n")
               
            if progress is not None:
                progress.update(n % 10000)
            infiles.close()
            # cutadapt 未提供接口可以记录
            # 追加输出到文件中，每个进程一个文件
            fname = f'{self.outdir}/_multithread_temp/barcodes_{pid}.txt'
            with open(fname, 'a') as f:
                f.writelines(results)

            return n,cbmatched,perfect_match_count

    with make_runner(input_paths, cores=cores) as runner:
        outfiles = NullOutputFiles()
        pipeline = SimplePipeline()
        stats = runner.run(pipeline, Progress(every=1), outfiles=outfiles)

    import subprocess
    final_path = f"{outdir}/spatial_umis.csv"

    # 拼接内容文件（按顺序追加）
    subprocess.run(f"cat {temp_dir}/barcodes_*.txt > {final_path} && rm -rf {temp_dir}", shell=True, check=True)

    results = pl.read_csv(f"{outdir}/spatial_umis.csv",
                          new_columns=["Cell_Barcode", "UMI", "Spatial_Barcode", "Read_Count"],has_header=False)
    results = results.group_by(["Cell_Barcode", "UMI", "Spatial_Barcode"]).agg(
        pl.col("Read_Count").sum()
    )
    results.write_csv(final_path)
    os.system(f"gzip {final_path}") ###压缩

    print(f"Total reads processed: {stats.n:,}, Matched reads: {stats.total_bp[0]:,}, Perfect matches: {stats.total_bp[1]:,}")
    return results,stats.n,stats.total_bp[0]


def calculate_statistics(df: pd.DataFrame, total_reads: int, matched_reads: int) -> Dict[str, float]:

    valid_reads = df["Read_Count"].sum()
    valid_umis = len(df)
    
    return {
        "Total Spatial Reads": total_reads,
        "Spatial Reads with Valid Cellbarcode": matched_reads,
        "Valid Spatial Reads": valid_reads,
        "Fraction of Valid Spatial Reads": round(valid_reads / total_reads, 4),
        "Valid Spatial UMIs": valid_umis,
        "Spatial Barcode Saturation": round(1 - (valid_umis / valid_reads), 4) if valid_reads else 0
    }


def stat_spatial_barcodes(r1fastq, r2fastq, 
                          oligochip, 
                          rnachemistry, mapparams, 
                          cbwhitelist, sbwhitelist, 
                          outdir, n_jobs) -> None:
    """主处理函数"""
    print("Starting spatial barcode analysis")
    os.makedirs(outdir, exist_ok=True)
    # 1. 准备白名单和位置信息
    pos_info = check_rnachemistry(rnachemistry, mapparams)
    cb_whitelist = get_whitelist(cbwhitelist, rnachemistry, mapparams)
    sb_whitelist = read_whitelist_files(sbwhitelist) if sbwhitelist else None
    cb_positions = list(zip(pos_info["cbstart"], pos_info["cblen"]))
    umi_pos = (pos_info["umistart"], pos_info["umilen"])
    # 2. 处理
    if n_jobs == 1:
        result_df, total_reads, total_matched = process_fastq_singlethread(
            r1fastq, r2fastq, oligochip,
            cb_positions=cb_positions,
            umi_pos=umi_pos,
            cb_whitelist= cb_whitelist, 
            sb_whitelist= sb_whitelist
        )
        spatial_umis_path = os.path.join(outdir, 'spatial_umis.csv')
        result_df.write_csv(spatial_umis_path)
    else:
        result_df, total_reads, total_matched = process_fastq_multithread(
            r1fastq, r2fastq, oligochip,
            cb_positions=cb_positions,
            umi_pos=umi_pos,
            cb_whitelist= cb_whitelist, 
            sb_whitelist= sb_whitelist,
            cores=n_jobs,
            outdir=outdir
        )

    # 3. 计算并保存统计信息
    summary = calculate_statistics(result_df, total_reads, total_matched)
    outstat = os.path.join(outdir, "sb_library_summary.temp.csv")
    with open(outstat, "wt") as outf:
        for k, v in summary.items():
            print(f"{k},{v}", file=outf)

    print(f"Statistics saved to {outstat}")
    del result_df

def parse_args():
    parser = argparse.ArgumentParser(description='Counting the spatial umi saturation')
    parser.add_argument('-r1', '--r1fastq', 
        metavar='FILE', 
        type=str,
        help='The R1 fastq file, contain cell barcode'
        )
    parser.add_argument('-r2', '--r2fastq', 
        metavar='FILE', 
        type=str,
        help='The R2 fastq file, contain spatial barcode and linkers'
        )
    parser.add_argument('-sb', '--sbwhitelist', 
        metavar='FILE', 
        type=str,
        help='The spatial barcode whitelist files'
        )
    parser.add_argument('-cb', '--cbwhitelist', 
        metavar='FILE', 
        type=str,
        help='The cell barcode whitelist files'
        )
    parser.add_argument('-rc', '--rnachemistry',
        metavar='STR',
        type=str,
        choices=['10X', 'leader_v1'],
        default='leader_v1',
        help='rnachemistry version: 10X or leader_v1, [default: leader_v1]'
        )
    parser.add_argument(
        '-oc', '--oligochip',
        metavar='STR',
        choices=["LD","GM"],
        help='spatial chip version, only can be "LD" or "GM". [default: LD].',
        default='LD'
        )
    parser.add_argument('-o', '--outdir', 
        metavar='PATH', 
        type=str,
        help='The output directory'
        )
    parser.add_argument(
        '-m', '--mapparams',
        metavar='STR',
        help='Additional STAR mapping parameters. must be provide while setting rnachemistry to other.'
        )
    parser.add_argument('-j', '--n_jobs', 
        metavar='INT', 
        type=int,
        help='Number of parallel workers, default: 4',
        default=4
        )
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_args()

    stat_spatial_barcodes(args.r1fastq, args.r2fastq, 
                          args.oligochip, 
                          args.rnachemistry, args.mapparams, 
                          args.cbwhitelist, args.sbwhitelist, 
                          args.outdir, args.n_jobs)