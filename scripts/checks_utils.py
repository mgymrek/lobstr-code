from subprocess import Popen, PIPE, STDOUT

import getopt
import os
import sys
###########################


def ExecuteCmd(cmd, debug):
    """ Execute a shell command and get stdout """
    if debug: 
        sys.stderr.write("CMD: " +cmd + "\n")
        return "debug"
    else:
        p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, \
                      stderr=STDOUT, close_fds=True)
        output = p.stdout.read()
        return output.strip()

def GetLineCount(filename, debug, zipped=False,):
    if zipped:
        cmd = "zcat %s | grep -v Period | grep -v version |  wc -l"%filename
        return ExecuteCmd(cmd, debug)
    else:
        cmd = "cat %s | grep -v Period | grep -v version |  wc -l"%filename
        return ExecuteCmd(cmd, debug)

def GetListFromStdout(stdout_string):
    """
    Convert column of items to a python list
    """
    return [item.strip() for item in stdout_string.split("\n")]
