import os
import sys
import json
import time
import logging
import sys
import io
import shutil
import base64
import subprocess
from datetime import datetime
from pathlib import Path
from datetime import timedelta
from functools import wraps
from typing import Union, List
from space_sketcher.__init__ import __root_dir__

# 配置高性能日志系统
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler(sys.stdout)]
)
logger = logging.getLogger(__name__)


def add_log(func):
    """
    logging start and done.
    """
    logFormatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    module = func.__module__
    name = func.__name__
    logger_name = f"{module}.{name}"
    logger = logging.getLogger(logger_name)
    logger.setLevel(logging.INFO)

    consoleHandler = logging.StreamHandler(sys.stderr) ##标准输出
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    @wraps(func) #在装饰器前加@wraps(func)能帮助保留原有函数的名称和文档字符串
    def wrapper(*args, **kwargs):
        start_time = time.time()  # 记录开始时间
        logger.info(f"{name} start...")
        result = func(*args, **kwargs) ###执行函数
        end_time = time.time()  # 记录结束时间

        run_time = str(round((end_time - start_time) / 60, 2))
        logger.info(f"{name} done. time used: {run_time}")
        return result 
    return wrapper

def csv2dict(filename: Path | str, sep=","):
    import pandas as pd
    
    df = pd.read_csv(filename, sep=sep, header=None)
    return dict(zip(df.iloc[:, 0], df.iloc[:, 1]))

# Extracted common path construction logic into a separate function
def get_common_path_part():
    # print(__root_dir__)
    return '/'.join(str(__root_dir__).split('/')[0:-1])

# Use os.makedirs with exist_ok=True instead of os.system
# Construct target path using os.path.join for better portability
# Handle exceptions and provide meaningful error messages
def str_mkdir(arg):
    """
    Exceptions:
    - If the directory already exists, no action is taken.
    - If there is no permission to create the directory, raises a PermissionError.
    - If any other OS error occurs, it is raised as-is.
    """
    try:
        os.makedirs(arg)  # Try creating the directory
    except FileExistsError:
        pass  # If the directory already exists, take no action
    except PermissionError:
        raise PermissionError("Permission denied to create the directory")  # If no permission to create, raise an exception
    except Exception as e:
        raise (f"Failed to create directory {arg}: {e}")  

# Use string formatting instead of concatenation for readability
# Handle exceptions and provide meaningful error messages    
def change_path():
    common_path = get_common_path_part()
    try:
        os.environ['PATH'] += f":{common_path}/bin" 
        os.environ['LD_LIBRARY_PATH'] = f"{common_path}/lib"
    except Exception as e:
        raise (f"Failed to update environment variables: {e}")  

# Construct bin path using os.path.join for better portability
def bin_path():
    # return os.path.join(get_common_path_part(), '.venv/bin')  
    return "/data03/lead/userdata/huanghuichang/Software/miniconda3/bin/"
    
def rm_temp(*args):
    """
    Remove specified files or directories, including their contents if they are directories.
    
    Parameters:
    *args (str): A variable-length argument list containing the file/directory paths to be removed.

    Note:
    - Symbolic links pointing to directories are skipped without attempting removal.
    - Any encountered exceptions during the removal process are caught and printed to stdout.
    """
    for filename in args:
        try:
            if os.path.exists(filename):
                if os.path.isdir(filename):
                    if os.path.islink(filename) or os.path.realpath(filename) != filename:
                        print(f"Skipped symbolic link: {filename}")
                        continue
                    shutil.rmtree(filename)
                else:
                    os.remove(filename)
            else:
                pass
        except Exception as e:
            print(f"Error removing {filename}: {e}")

def get_formatted_time():
    current_time = datetime.now()
    formatted_time = current_time.strftime('%Y-%m-%d %H:%M:%S')
    return formatted_time


def start_print_cmd(
    command: Union[str, List[str]],  # 支持字符串或列表形式的命令
    name: str,                      # 日志标识（如模块名）
    log_dir: str,                   # 日志目录
    ): 
    logger = logging.getLogger()
    logger.info(command)
    # 判断命令是否需要执行
    if not os.path.exists(os.path.join(log_dir, f".{name}.done")) or not os.path.exists(os.path.join(log_dir, f".{name}.skipped")):
        # 执行命令并实时输出
        try:
            subprocess.check_call(command, shell=True)
        except subprocess.CalledProcessError as e:
            logger.error("[FAILED] Exit code: %d\nCommand: %s", 
                        e.returncode, command)
            raise
        except Exception as e:
            logger.error("[UNEXPECTED ERROR] %s", str(e))
            raise
        # # 创建完成标记文件
        # with open(os.path.join(log_dir, f".{name}.done"), 'w') as f:
        #     f.write("done")
    else:
        logger.info(f".{name}.done or .{name}.skipped already exists, skip running {name}")


# def setup_logging(name, log_dir):
#     today = time.strftime('%Y%m%d', time.localtime(time.time()))
#     logfile = f'{log_dir}/log/{today}.txt'
#     logger = logging.getLogger(name)
    
#     if not logger.handlers:
#         logger.setLevel(logging.INFO)
#         formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s \n%(message)s')
        
#         file_handler = logging.FileHandler(logfile, encoding="utf8")
#         file_handler.setLevel(logging.DEBUG)
#         file_handler.setFormatter(formatter)
        
#         console_handler = logging.StreamHandler(sys.stdout)
#         console_handler.setLevel(logging.ERROR)
#         console_handler.setFormatter(formatter)
        
#         logger.addHandler(file_handler)
#         logger.addHandler(console_handler)
    
#     return logger

# def logging_call(popenargs, name, log_dir):
#     logger = setup_logging(name, log_dir)

#     try:
#         output = subprocess.check_output(popenargs, shell=True, stderr=subprocess.STDOUT, universal_newlines=True)
#         logger.info('%s', output)
#     except subprocess.CalledProcessError as e:
#         logger.error('Command failed with exit code %d', e.returncode)
#         logger.error('%s', e.output)


def judgeFilexits(*args):
    # Flatten args and filter out empty strings
    files_to_check = [file for arg in args for file in arg.split(',') if file]
    # Check for empty file names
    if any(not file for file in files_to_check):
        print("Error: Received empty file name(s).")
        return
    # Use a set comprehension to check file existence and collect missing files
    missing_files = {file for file in files_to_check if not os.path.exists(file)}
    # If any files are missing, print error messages and raise a custom exception
    if missing_files:
        error_msgs = [" ------------------------------------------------"]
        for file in missing_files:
            error_msgs.append("Error: Cannot find input file or directory '{}'".format(file))
        error_msgs.append(" ------------------------------------------------")
        print("\n".join(error_msgs), end="\n\n")
        raise FileNotFoundError("One or more input files do not exist.")


def gunzip(in_file):
    """
    命令行 gunzip 命令
    """
    import subprocess
    import tempfile

    try:
        with tempfile.NamedTemporaryFile(delete=False, suffix=".txt") as temp_file:
            subprocess.check_call(f"gunzip -c {in_file} > {temp_file.name}", shell=True)
            return temp_file.name
    except Exception as e:
        raise ValueError(f"解压缩文件出错: {e}")
