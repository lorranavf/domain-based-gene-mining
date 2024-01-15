import importlib
import subprocess
def check_libraries():

    libraries = ['biolib', 'Bio', 'pandas', 'seaborn']

    missing_libraries = []
    for library in libraries:
        try:
            __import__(library, fromlist=['']).__name__
        except ModuleNotFoundError:
            missing_libraries.append(library)

    if missing_libraries:
        print('The following libraries need to be installed:')
        for library in missing_libraries:
            print(library)
    else:
        print('All the required libraries are installed.')


def check_programs():

    programs = ['docker', 'mamba', 'makeblastdb', 'blastp', 'hmmscan', 'cath-resolve-hits', 'signalp6',
                'pepstats', 'deeploc2', 'mafft', 'CIAlign', 'iqtree', 'orthofinder', 'ete3']

    missing_programs = []
    
    for program in programs:
        try:
            # Execute the shell command to check if the program is present
            subprocess.check_output(['which', program], stderr=subprocess.DEVNULL)

        except subprocess.CalledProcessError:
            missing_programs.append(program)

    if missing_programs:
        print('The following programs are missing in the environment:')
        for program in missing_programs:
            print(program)
    else:
        print('All the required programs are present in the environment.')

check_libraries()
check_programs()
