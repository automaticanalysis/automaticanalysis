# Note: Assume Python 2.7, since that will have been installed for facemasking requirements.
#
# Using python (instead of matlab) because using system() in matlab gives enormous amount of
# command window output, and there is some problem with the status code from system()

import sys
import subprocess
import platform

def main():
    do_install()

def do_install():
    plt = platform.system()
    if plt == "Windows":
        print("Your system is Windows")
        sys.stdout.flush()
        # do x y z
    elif plt == "Linux":
        print("Your system is Linux")
        sys.stdout.flush()
        subprocess.check_call('source $GITHUB_WORKSPACE/.github/workflows/tools_install.sh', shell=True, executable="/bin/bash")

    elif plt == "Darwin":
        print("Your system is MacOS")
        sys.stdout.flush()
        # do x y z
    else:
        print("Unidentified system:")
        print(plt)
        sys.stdout.flush()

if __name__=="__main__":
    main()