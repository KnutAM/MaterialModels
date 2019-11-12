import os
import sys


def main(input_args):

    ftyps = ['RF', 'dRdX']  # The two function types that should exist
        
    (umod, ftyp, nbak) = get_flist(ftyps)

    
    for m,f,n in zip(umod, ftyp, nbak):
        write_module(m,f,n)
    
    
def get_flist(ftyps):
    dirs = os.listdir('build/')
    
    ftypes = []
    models = []
    nbacks = []
    for file in dirs:
        tmp = file.split('.')
        if tmp[-1] == 'f90':
            base = tmp[0].split('_')
            if base[0] in ftyps:
                ftypes.append(base[0])
                models.append('_'.join(base[1:-1]))
                nbacks.append(int(base[-1]))
            else:
                print('Unknown function type ' + base[0])
    
    # Organize
    umodels = list(set(models))
    ftyp = []
    nbak = []
    for um in umodels:
        ftyp.append([])
        nbak.append([])
        for m,f,n in zip(models, ftypes, nbacks):
            if m==um:
                ftyp[-1].append(f)
                nbak[-1].append(n)
                
    return umodels, ftyp, nbak
    

def write_module(model, ftypes, nback):
    module_name = 'build_output/' + model + '_acegen_mod.f90'
    fid = open(module_name, 'w')
    fid.write('! Module for ' + model + '\n')
    fid.write('! Functions generated with AceGen and collected by make_acegen_mods.py' + '\n')
    fid.write('module acegen_mod' + '\n')
    fid.write('use smsutility' + '\n')
    fid.write('implicit none' + '\n\n')
    fid.write('private' + '\n\n')
    for f,n in zip(ftypes, nback):
        fid.write('public ' + f + str(n) + '\n')
    
    fid.write('\n')    
    fid.write('contains' + '\n\n')
    
    for f,n in zip(ftypes, nback):
        fname = 'build/' + f + '_' + model + '_' + str(n) + '.f90'
        fid.write(clean_filestring(fname))
        
    fid.write('end module acegen_mod')
    fid.close()
    

def clean_filestring(file):
    # Remove comments above subroutine
    # Remove the v argument from the subroutine
    # Remove the use SMSUtility line
    # Add " SUBROUTINE" to the last "END" statement (Needed?) (not implemented...)
    
    with open(file, 'r') as fid:
        for line in fid:
            if len(line.split('!')[0])>1:
                break # First non-empty, non-commented line found
        line = '('.join(line.split('(v,'))
        fstr = line
        line = fid.readline()   # Discard the next line (use SMSUtility)
        
        for line in fid:
            fstr = fstr + line
        
        fstr = fstr[:-1] + ' subroutine' + '\n\n'    # Add end subroutine statement
        
    return fstr
        

if __name__ == '__main__':              
    main(sys.argv)                              # run the main function
