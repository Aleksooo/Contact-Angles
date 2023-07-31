import numpy as np
import argparse
import os
from src.utils_py.io.gro import read_gro, write_gro
from src.utils_py.geom.Box import Box
from src.utils_py.assembler.build import build_system

parser = argparse.ArgumentParser()
parser.add_argument('H', type=float, help='Slot H')
parser.add_argument('phi', type=float, help='Amoung of water')
parser.add_argument('build', type=int, help='Build system flag')
parser.add_argument('path', type=str, help='path on server')
parser.add_argument('subf', type=str, help='folder with substrate')

args = parser.parse_args()
structure = read_gro(f'{args.subf}/cal.gro')

WIDTH_X, WIDTH_Y, HEIGHT = structure.box
# WIDTH_X = 9.23520
# WIDTH_Y = 13.49130
# # H = WIDTH_Y / 3 # 5.34

H = args.H

offset = 2
delta_h = 0.2
frac = args.phi

# Конфигурация половины
folder = f'cal_{H:.1f}_{frac:.1f}_{HEIGHT:.1f}'
dir = 'ff/gromos'

''' ПАРАМЕТРЫ НЕ ТРОГАТЬ!!! '''
insertion_limit = int(1e5)
rotation_limit = 1000
package = 0.3
distance = {'min': 0.08**2, 'opt': 0.12**2}

system_size = np.array([WIDTH_X, WIDTH_Y, HEIGHT + H])
# structure = read_gro(f'{dir}/gro/cal_fixed.gro')
points = list(structure.get_XYZ())
structure.box = system_size

names = ['decane', 'spce']
density = [3, 33] # nm-3

box_left = Box(
    center=np.array([WIDTH_X/2, WIDTH_Y*(1-frac)/2, HEIGHT + H/2]),
    borders=np.array([WIDTH_X, WIDTH_Y*(1-frac), H])
)

box_left_delta = Box(
    center=np.array([WIDTH_X/2, WIDTH_Y*(1-frac)/2, HEIGHT + H/2]),
    borders=np.array([WIDTH_X, WIDTH_Y*(1-frac), H - 2 * delta_h])
)

box_right = Box(
    center=np.array([WIDTH_X/2, WIDTH_Y*(1-frac/2), HEIGHT + H/2]),
    borders=np.array([WIDTH_X, WIDTH_Y*frac, H])
)

box_right_delta = Box(
    center=np.array([WIDTH_X/2, WIDTH_Y*(1-frac/2), HEIGHT + H/2]),
    borders=np.array([WIDTH_X, WIDTH_Y*frac, H - 2 * delta_h])
)


insert_shapes = [box_left_delta, box_right_delta]
shapes = [box_left, box_right]
numbers = list(np.round(np.array([shapes[i].get_volume() * density[i] for i in range(len(names))])).astype(int))

if args.build:
    structure = build_system(
        dir, structure, names, numbers, insert_shapes, points,
        insertion_limit=insertion_limit,
        rotation_limit=rotation_limit,
        package=package,
        min_dist2=distance['min']
    )


# Запись .gro и system.itp
mol_names = ''.join(map(lambda x: 'w' if x == 'spce' else x[0], names))
mol_nums = '_'.join(map(str, numbers))
filename = 'cal_'+ mol_names + '_' + mol_nums

os.system(f'mkdir na_diff/{folder}')

with open(f'na_diff/{folder}/system.itp', 'w') as f:
    for name in names:
        f.write(f'#include "../{name}.itp"\n')
    f.write('#include "cal.itp"\n')

    f.write(f'\n[ system ]\n{filename}\n')
    f.write('\n[ molecules ]\n; molecule name\tnr.\n')
    f.write('#include "calmol.itp"\n')
    for i, name in enumerate(names):
        f.write(f'{name}\t{numbers[i]}\n')

print('Writing .gro files.')

if args.build:
    with open(f'na_diff/{folder}/{filename}.gro', 'w') as f:
        f.write(write_gro(structure))

    # writing backup
    with open(f'na_diff/{folder}/#{filename}.gro#', 'w') as f:
        f.write(write_gro(structure))

    # Перемешивание
    print('Mixing system')
    os.system(f'./mixer -f na_diff/{folder}/{filename}.gro -o na_diff/{folder}/{filename}.gro -mn2 {distance["min"]} -opt2 {distance["opt"]}')

# Создание скрипта для запуска
for sc in [1]:
    with open(f'na_diff/{folder}/run_x{sc}.sh', 'w') as f:
        f.write(f'#!/bin/bash\n#SBATCH --job-name=gromacs_{filename}\n#SBATCH --partition=RT\n#SBATCH --nodes=3\n#SBATCH --ntasks-per-node=16\n\n')
        f.write(f'cp ../gromos_x{sc}.top gromos_x{sc}.top\n')
        f.write(f'srun -n 1 gmx_mpi grompp -f ../nvt_cal_steep.mdp -c {filename}_init.gro -p gromos_x{sc}.top -o {filename} -maxwarn 10\n')
        f.write(f'srun -n 16 gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm {filename}\n')
        f.write(f'srun -n 1 gmx_mpi grompp -f ../nvt_cal_short.mdp -c {filename}.gro -p gromos_x{sc}.top -o {filename} -maxwarn 10\n')
        f.write(f'srun -n 48 gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm {filename} -dlb yes\n')
        f.write(f'srun -n 1 gmx_mpi grompp -f ../nvt_cal_run.mdp -c {filename}.gro -p gromos_x{sc}.top -o {filename} -maxwarn 10\n')
        f.write(f'srun -n 48 gmx_mpi mdrun -s -o -x -c -e -g -v -deffnm {filename} -dlb yes\n')
        f.write(f'rm ./#*#')

# Отправка файлов на сервер при перезапуске
print(f'Sending files to server to the folder: "{folder}"')

os.system(f"ssh mipt-nd 'mkdir {args.path}'")
os.system(f"ssh mipt-nd 'mkdir {args.path}/{folder}'")

files = ['spce.itp', 'decane.itp', 'cal.itp']
for file in files:
    os.system(f'scp {dir}/itp/{file} mipt-nd:{args.path}/{file}')

os.system(f'scp {dir}/itp/cal.itp mipt-nd:{args.path}/{folder}/cal.itp')
os.system(f'scp {args.subf}/calmol.itp mipt-nd:{args.path}/{folder}/calmol.itp')

files = ['nvt_cal_steep.mdp', 'nvt_cal_short.mdp', 'nvt_cal_run.mdp']
for file in files:
    os.system(f'scp {dir}/mdp/{file} mipt-nd:{args.path}/{file}')

files = [filename+'.gro', 'system.itp', 'run_x1.sh']
for file in files:
    os.system(f'scp na_diff/{folder}/{file} mipt-nd:{args.path}/{folder}/{file}')
os.system(f'scp na_diff/{folder}/{filename}.gro mipt-nd:{args.path}/{folder}/{filename}_init.gro')
