import os
import sys
import math
import typer
from typing_extensions import Annotated
from typing import Optional
from subprocess import check_call
from space_sketcher.tools.utils import judgeFilexits, change_path
from space_sketcher.__init__ import __root_dir__

def count_chromosomes(genome_file):
    chromosomes = set()
    with open(genome_file, "r") as fasta_file:
        for line in fasta_file:
            if line.startswith(">"):
                chromosome = line.strip().lstrip(">")
                chromosomes.add(chromosome)
    return len(chromosomes)

def star_index(fasta,gtf,genomeDir,star_program,limitram,threads):
    if not os.path.exists(genomeDir):
        os.system('mkdir -p %s'%genomeDir)
    genome_size = os.path.getsize(fasta)
    SAindexNbases = min(int(math.log2(genome_size)/2 - 1),14)
    n_entries = count_chromosomes(fasta)
    chr_bins = min(18, int(math.log2(max(int(genome_size / n_entries),100))))
    
    star_cmd = [
        star_program,
        '--runMode genomeGenerate',
        f'--runThreadN {threads}',
        f'--genomeDir {genomeDir}',
        f'--genomeFastaFiles {fasta}',
        f'--sjdbGTFfile {gtf}',
        '--sjdbOverhang 99',
        f'--limitGenomeGenerateRAM {limitram}',
        f'--genomeSAindexNbases {SAindexNbases}',
        f'--genomeChrBinNbits {chr_bins}'
    ]
    star_cmd_str = ' '.join(star_cmd)

    print('STAR verison: 2.7.11b')
    print('runMode: genomeGenerate')
    print('runThreadN: %s'%threads)
    print('limitGenomeGenerateRAM: %s'%limitram)
    print('genomeSAindexNbases: %s'%SAindexNbases)
    print('genomeChrBinNbits: %s'%chr_bins)
    print('genomeDir: %s'%genomeDir)
    print('fasta: %s'%fasta)
    print('gtf: %s'%gtf)

    sys.stdout.flush()
    check_call(star_cmd_str,shell=True)


class Ref:
    def __init__(self, args):
        self.ingtf = os.path.abspath(args.ingtf)
        self.fasta = os.path.abspath(args.fasta)
        self.species = args.species
        self.genomeDir = os.path.abspath(args.genomeDir)
        self.limitram = args.limitram
        self.threads = args.threads

    def run(self):
        change_path()
        judgeFilexits(self.ingtf,self.fasta)
        self.genomeDir = os.path.abspath(self.genomeDir)
        starbin = os.path.join(__root_dir__, "software/STAR")
        if not self.noindex:
            star_index(self.fasta,
                       self.ingtf,
                       self.genomeDir,
                       starbin, 
                       self.limitram, 
                       self.threads)
        print("\033[0;32;40mAnalysis Complete\033[0m")


def mkref_app(
    # 必需参数
    fasta: Annotated[str, typer.Option("--fasta", "-f", help="Path to the genome file in FASTA format", prompt=True, show_default=False)],
    ingtf: Annotated[Optional[str], typer.Option("--ingtf", "-i", help="Path to the genome annotation file in GTF format", prompt=True, show_default=False)],
    # 可选参数
    species: Annotated[str, typer.Option("--species", "-s", help="Species name ")] = "undefined",
    genomeDir: Annotated[str, typer.Option("--genomeDir", "-g", help="Path to store database files ")] = os.getcwd(),
    limitram: Annotated[int, typer.Option("--limitram", "-l", help="Maximum RAM (bytes) for indexing ")] = 125000000000,
    threads: Annotated[int, typer.Option("--threads", "-t", help="Number of threads ")] = 4
):
    """
    Build RNA reference database.\n
    
    Example:\n
        space-sketcher rna mkref --fasta /path/to/fasta --ingtf /path/to/ingtf --species human --genomeDir /path/to/genomeDir
    """
    # 将参数转换为类似argparse的Namespace对象
    class Args:
        pass
    args = Args()
    for k, v in locals().items():
        if k != 'args':
            setattr(args, k, v)
    
    # 执行处理流程
    processor = Ref(args)
    processor.run()

# 导出函数
__all__ = ["mkref_app"]

