import textwrap
from colorama import init, Fore, Style

# 初始化 colorama
init()

def sum_help():
    return textwrap.dedent(f'''\
        {Fore.CYAN}Spatial Transcriptomic Analysis Workflow Suite:{Style.RESET_ALL}
        --------------------------------
        {Fore.GREEN}space_sketcher rna{Style.RESET_ALL}      {Fore.YELLOW}Spatial RNA Analysis Workflow{Style.RESET_ALL}
        This suite provides a comprehensive set of tools for spatial transcriptomic analysis, including RNA sequencing workflows.
        ''')

def help_text(category, prefix):
    if category == 'rna':
        return textwrap.dedent(f'''\
            {Fore.CYAN}Spatial RNA Analysis Workflow:{Style.RESET_ALL}

            {Fore.GREEN}Function:{Style.RESET_ALL}
            {Fore.GREEN}{prefix} run{Style.RESET_ALL}
                     Utilize spatial cDNA and spatial oligo library sequencing data for quality control,   
                     alignment, and functional region annotation analysis. Detect droplets containing 
                     cells to generate a filtered gene expression matrix. 
                     Apply the filtered gene 
                     expression matrix for cell filtering, dimensionality  
                     reduction, clustering, and annotation analysis.                                       

            {Fore.GREEN}{prefix} multi{Style.RESET_ALL}
                     Generate shell scripts to execute the 'run' command on multiple samples.              

            {Fore.GREEN}{prefix} mkref{Style.RESET_ALL}
                     Construct reference database.
            ''')
    # elif category == '':
    # more and more function to be added.
    else:
        return f"Unknown category: {category}"



# def sum_help():
#     return textwrap.dedent(f'''\
#         Spatial Transcriptomic Analysis Workflow Suite:
#         --------------------------------
#         space_sketcher rna      Spatial RNA Analysis Workflow
#         This suite provides a comprehensive set of tools for spatial transcriptomic analysis, including RNA sequencing workflows.
#         ''')

# def help_text(category, prefix):
#     if category == 'rna':
#         return textwrap.dedent(f'''\
#             Spatial RNA Analysis Workflow:

#             Function:
#             {prefix} run
#                      Utilize spatial cDNA and spatial oligo library sequencing data for quality control,   
#                      alignment, and functional region annotation analysis. Detect droplets containing 
#                      cells to generate a filtered gene expression matrix. 
#                      Apply the filtered gene 
#                      expression matrix for cell filtering, dimensionality  
#                      reduction, clustering, and annotation analysis.                                       

#             {prefix} multi
#                      Generate shell scripts to execute the 'run' command on multiple samples.              

#             {prefix} mkref
#                      Construct reference database.
#             ''')
#     # elif category == '':
#     # more and more function to be added.
#     else:
#         return f"Unknown category: {category}"