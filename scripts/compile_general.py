from __future__ import print_function
import sys
import os
import shutil
import platform
from pathlib import Path

if platform.system() == 'Linux':
    linux_compiler = 'ifort'
    dll_suffix = '.so'
    dll_location = ''
    cmake_setup = 'FC=' + linux_compiler + ' cmake ../src'
elif platform.system() == 'Windows':
    dll_suffix = '.dll'
    dll_location = 'Release/'
    cmake_setup = 'cmake -G "Visual Studio 12 2013 Win64" ../src'
elif platform.system() == 'Darwin':
    print('MAC not supported by this script')
else:
    print('Unknown platform (os), not supported by this script')

all_model_categories = ['GenSmallStrain', 'GenFiniteStrain']

output_dir = Path('../compiled/').absolute()
script_dir = Path(os.path.realpath(__file__)).parent

if not output_dir.exists():
    os.mkdir(output_dir)


def main(argv):
    # Input argument 1
    #   Which model category to compile (default is all_model_categories)
    # Input arguments 2-
    #   Sequence of indicies (row nr - 1) of models to compile, see cmake_compile_models.txt or abaqus_compile_files.txt
    #   if the former doesn't exist. The default is to compile all models in category

    if len(argv) > 1:
        model_categories = [argv[1]]
    else:
        model_categories = all_model_categories

    if len(argv) > 2:
        line_numbers = []
        for nr in argv[2:]:
            line_numbers.append(int(nr))
    else:
        line_numbers = None

    for model_cat in model_categories:
        os.chdir(script_dir)
        build_model_cat(model_cat, line_numbers)


def build_model_cat(model_cat, line_numbers=None):
    # Tasks
    # For each row in model_cat/src/abaqus_compile_files.txt
    #   Call routine to build specific model
    # If line_numbers is not None, only call build routine at those particular lines
    model_cat_folder = Path('../models/' + model_cat)

    models, cmake_set_vars = get_models(model_cat_folder, line_numbers)

    os.chdir(model_cat_folder)

    build_dir = 'build/'

    if os.path.exists(build_dir):
        if input('Do you want to clean build? (yes/no)= ') == 'yes':
            shutil.rmtree(build_dir)
            os.mkdir(build_dir)
    else:
        os.mkdir(build_dir)

    os.chdir(build_dir)

    os.system(cmake_setup)
    for mod, var in zip(models, cmake_set_vars):
        print('cmake ' + var + ' .')
        os.system('cmake ' + var + ' .')
        os.system('cmake --build . --config RELEASE')
        shutil.copy(dll_location + mod + dll_suffix, output_dir)


def get_models(model_cat_folder, line_numbers):
    custom_cmake_specs = model_cat_folder / 'src' / 'cmake_compile_models.txt'

    cmake_var_set = []
    models = []
    if custom_cmake_specs.exists():
        specs = custom_cmake_specs.read_text().splitlines()
        if line_numbers is None:
            for spec in specs:
                models.append(spec.split(':')[0])
                cmake_var_set.append(spec.split(':')[1])
        else:
            for nr in line_numbers:
                if nr >= len(specs):
                    print(model_cat_folder.name + ' has only ' + str(len(specs)) + ' specific models')
                    print('skipping nr ' + str(nr + 1))
                else:
                    models.append(specs[nr].split(':')[0])
                    cmake_var_set.append(specs[nr].split(':')[1])
    else:
        lines = (model_cat_folder / 'src' / 'abaqus_compile_files.txt').read_text().splitlines()
        if line_numbers is None:
            for line in lines:
                models.append(line.split(',')[0])
                cmake_var_set.append('-D model=' + models[-1])
        else:
            for nr in line_numbers:
                if nr >= len(lines):
                    print(model_cat_folder.name + ' has only ' + str(len(lines)) + ' specific models')
                    print('skipping nr ' + str(nr+1))
                else:
                    models.append(lines[nr].split(',')[0])
                    cmake_var_set.append('-D model=' + models[-1])

    return models, cmake_var_set


if __name__ == '__main__':
    main(sys.argv)  # run the main function
