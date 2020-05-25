import sys
import re
from pathlib import Path
import os



def main(argv):
    print(os.getcwd())
    module_name = "acegen_crystal_plasticity_mod"
    subroutine_names = ("model_size", "elastic", "residual", "jacobian", "plastic_output")
    if len(argv) > 1:
        name = argv[1]
    else:
        name = "tmp"

    with open(name + ".f90", "w") as fid:
        write_header(fid, module_name)
        for name in subroutine_names:
            sub = fix_subroutine(name)
            fid.write(sub)
        fid.write("\n" + "end module " + module_name)


def fix_subroutine(name):
    subroutine = Path(name+'.f90').read_text()
    subroutine = re.sub("SUBROUTINE " + name + "\(v,", "SUBROUTINE " + name + "(", subroutine)
    subroutine = re.sub("END$", "END SUBROUTINE", subroutine)

    return subroutine


def write_header(fid, module_name):
    with open('tmpinfo.txt') as inp_fid:
        for line in inp_fid:
            fid.write('!' + line)

    fid.write("\n" + "module " + module_name + "\n")
    fid.write("implicit none" + "\n")
    fid.write("contains" + "\n")


if __name__ == '__main__':
    main(sys.argv)
