from shutil import copyfile
import os

work_dir = os.getcwd()

path_parent = os.path.dirname(os.getcwd())
os.chdir(path_parent) 


if not os.path.exists(f'{work_dir}/input'):
    os.makedirs(f'{work_dir}/input') 
if not os.path.exists(f'{work_dir}/output'):
    os.makedirs(f'{work_dir}/output')


copy_from = 'step1_CreateNetworks/output'

copy_to = f'{work_dir}/input'
to_copy = ['COEX', 'GI', 'MI_MetAbundWeight', 'PPI_ppTW']

cell_conds = ['ND_18', 'WT_18', 'WT_32', 'ND_32', 'WT_8', 'ND_8', 'ND_25', 'WT_25', 'WT_37', 'ND_37']
for cc in cell_conds:
    save_dir = f'{copy_to}/{cc}'
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    for tc in to_copy:
        copyfile(f'{copy_from}/{cc}/{tc}_{cc}.edgelist', f'{save_dir}/{tc}_{cc}.edgelist')
    copyfile(f'{copy_from}/{cc}/E_{cc}_normalized.csv', f'{save_dir}/E_{cc}_normalized.csv')

os.chdir(work_dir) 
