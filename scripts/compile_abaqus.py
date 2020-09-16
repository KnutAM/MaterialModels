import os
import sys
import pathlib as pl
import shutil
import distutils.dir_util as dutil

all_model_categories = ['GenSmallStrain', 'GenFiniteStrain']

def main(argv):
    # Tasks, for each folder (model category):
    # For each line in "abaqus_compile_files.txt"
    # - Create a file umat.for or umat.f for windows or linux, respectively, with all content of included files
    # - Call abaqus make routine, i.e., abaqus make library="umat.[for/f]" (Note that previous builds should be removed)
    # - Rename output according to info in "abaqus_compile_files.txt" and move to the compiled_abaqus folder (if not
    #   existent, create this folder)

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
        build_model_cat(model_cat, line_numbers)


def build_model_cat(model_cat, line_numbers=None):
    # Tasks
    # For each row in model_cat/src/abaqus_compile_files.txt
    #   Call routine to build specific model
    # If line_numbers is not None, only call build routine at those particular lines
    model_cat_folder = pl.Path('../models/' + model_cat)
    lines = (model_cat_folder / 'src' / 'abaqus_compile_files.txt').read_text().splitlines()
    if line_numbers is None:
        for line in lines:
            build_model(model_cat, line)
    else:
        for nr in line_numbers:
            build_model(model_cat, lines[nr])


def build_model(model_cat, model_info):
    # Given the model_cat (giving the base folder containing the src folder) and
    # the model_info (the line in abaqus_compile_files.txt) containing the model name followed by each used file in
    # correct order (relative path from the src folder)
    # Combine all files into one umat and build it
    model_name = model_info.split(',')[0]
    model_files = model_info.split(',')[1:]

    print('Building model ' + model_cat + ' (' + model_name + ')')

    umat = create_umat()
    file_contents = ''
    for model_file in model_files:
        file = pl.Path('../models/' + model_cat) / 'src' / model_file.strip('"')
        file_contents = file_contents + file.read_text() + '\n'

    # dec_freeform (if written, is not needed within the file). Remove it, and add it in the beginning of the whole file
    dec_freeform = '!DEC$ FREEFORM'
    ind = file_contents.find(dec_freeform)
    file_contents = file_contents[:ind] + '\n' + file_contents[(ind+len(dec_freeform)):]
    file_contents = dec_freeform + '\n' + file_contents

    # Write to the umat.for / umat.f file
    umat.write_text(file_contents)

    dutil._path_created = {}  # Need to clear cache after deleting the old tree using shutil.rmtree
    dutil.copy_tree('../models/umat_utils', str(umat.parent))
    
    os.chdir(umat.parent)
    os.system('abaqus make library=' + str(umat.name))

    os.chdir('..')
    print(os.getcwd())

    compiled_dir = pl.Path('../compiled_abaqus')
    if not compiled_dir.exists():
        os.mkdir(compiled_dir)

    if str(umat).endswith('.for'):
        delete_files([model_name+'.dll', model_name+'.obj'], compiled_dir)
        os.rename(umat.parent / 'standardU.dll', compiled_dir / (model_name + '.dll'))
        os.rename(umat.parent / 'umat-std.obj', compiled_dir / (model_name + '.obj'))
    else:
        delete_files([model_name+'.so', model_name+'.o'], '../compiled_abaqus/')
        os.rename(umat.parent / 'libstandardU.so', compiled_dir / (model_name + '.so'))
        os.rename(umat.parent / 'umat-std.o', compiled_dir / (model_name + '.o'))
    
    # Move umat sources to a folder such that umat can be compiled together with other user 
    # subroutines if required
    src_dir = compiled_dir / model_name
    if src_dir.exists():
        shutil.rmtree(src_dir)
        
    os.rename(umat.parent, src_dir)
    

def create_umat():  # Umat actually not created, but path object to it is.
    dir = pl.Path('tmp_build_dir')
    if dir.exists():
        shutil.rmtree(str(dir))

    os.mkdir(str(dir))

    if sys.platform.startswith('win'):
        umat = dir / 'umat.for'
    elif sys.platform.startswith('linux'):
        umat = dir / 'umat.f'
    else:
        print('System not recognized')

    return umat


def delete_files(files, prepend=''):
    prepend = pl.Path(prepend)
    for file in files:
        if os.path.exists(prepend / file):
            os.remove(prepend / file)


if __name__ == '__main__':
    main(sys.argv)
