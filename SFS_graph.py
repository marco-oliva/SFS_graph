#!/usr/bin/env python3

import sys, time, argparse, subprocess, os, signal, shutil

Description = """
Rimerge
"""

#------------------------------------------------------------
def remove_file(file_path):
    os.remove(file_path)

def remove_dir(dir_path):
    shutil.rmtree(dir_path)

def move_dir_content(dir_source, dir_dest):
    files = os.listdir(dir_source)
    for f in files:
        shutil.move(dir_source + "/" + f, dir_dest)

#------------------------------------------------------------
# execute command: return True if everything OK, False otherwise
def execute_command(command, seconds):
    try:
        process = subprocess.Popen(command.split(), preexec_fn=os.setsid)
        process.wait(timeout=seconds)
    except subprocess.CalledProcessError:
        print("Error executing command line:")
        print("\t"+ command)
        return False
    except subprocess.TimeoutExpired:
        os.killpg(os.getpgid(process.pid), signal.SIGTERM)
        print("Command exceeded timeout:")
        print("\t"+ command)
        return False
    return True

