import numpy as np
import argparse
import os
from tqdm import tqdm
from src.utils_py.io.gro import read_gro, write_gro
from src.utils_py.geom.Roll import Roll
from src.utils_py.geom.AntiRoll import AntiRoll
from src.utils_py.assembler.build import build_system

parser = argparse.ArgumentParser()
parser.add_argument('H', type=float, help='Slot H')
parser.add_argument('phi', type=float, help='Amoung of water')
parser.add_argument('angle', type=int, help='angle of droplet')
parser.add_argument('subf', type=str, help='folder with substrate')
parser.add_argument('iter', type=int, help='folder with substrate')

args = parser.parse_args()
structure = read_gro(f'{args.subf}/cal.gro')

WIDTH_X, WIDTH_Y, HEIGHT = structure.box
H = args.H

offset = 2
delta_h = 0.2
frac = args.phi

# Конфигурация
folder = f'angles/{args.angle}'
dir = 'ff/trappe'

''' ПАРАМЕТРЫ НЕ ТРОГАТЬ!!! '''
insertion_limit = int(1e5)
rotation_limit = 1000
package = 0.3
distance = {'min': 0.08**2, 'opt': 0.12**2}

system_size = np.array([WIDTH_X, WIDTH_Y, HEIGHT + H])
print('Angle:', args.angle)
for iter in tqdm(range(args.iter)):
    structure = read_gro(f'{args.subf}/cal.gro')
    points = structure.atoms_xyz
    structure.box = system_size

    names = ['decane', 'tip4p']
    density = [3, 33]

    roll = Roll(
        center=np.array([WIDTH_X/2, WIDTH_Y/2, HEIGHT + H/2]),
        phi=0.5,
        theta=np.deg2rad(args.angle),
        borders=np.array([WIDTH_X, WIDTH_Y, H]),
    )

    anti_roll = AntiRoll(
        center=np.array([WIDTH_X/2, WIDTH_Y/2, HEIGHT + H/2]),
        phi=0.5,
        theta=np.deg2rad(args.angle),
        borders=np.array([WIDTH_X, WIDTH_Y, H]),
    )

    insert_shapes = [roll, anti_roll]
    shapes = [roll, anti_roll]
    numbers = list(np.round(np.array([shapes[i].get_volume() * density[i] for i in range(len(names))])).astype(int))

    structure = build_system(
        dir, structure, names, numbers, insert_shapes, points,
        insertion_limit=insertion_limit,
        rotation_limit=rotation_limit,
        package=package,
        min_dist2=distance['min']
    )

    # Запись .gro и system.itp
    filename = f'cal_{args.angle}_{iter+1}'

    os.system(f'mkdir {folder}')

    with open(f'{folder}/{filename}.gro', 'w') as f:
        f.write(write_gro(structure))

    # writing backup
    # with open(f'na_diff/{folder}/#{filename}.gro#', 'w') as f:
    #     f.write(write_gro(structure))

    # Перемешивание
    print('Mixing system')
    os.system(f'./mixer -f {folder}/{filename}.gro -o {folder}/{filename}.gro -mn2 {distance["min"]} -opt2 {distance["opt"]}')
