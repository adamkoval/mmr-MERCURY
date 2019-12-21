# Function file for simulation package

def read_envfile(envfile):
    """
    Function for reading envfile.
    """
    with open(envfile, 'r') as f:
        lines = [line for line in f.readlines() if line[0] is not '#']
        pyenv = lines[0].split()[2]
    return pyenv
