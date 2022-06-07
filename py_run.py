import os
import sys

if __name__ == "__main__":

    os.system("rm *.o")
    os.system("rm run_capelec")

    os.system("make")

    if len(sys.argv)>1 and sys.argv[1] == 'log':
        os.system("./run_capelec 1> logs/term_out.txt 2> logs/term_err.txt")
    else:
        os.system("./run_capelec")
