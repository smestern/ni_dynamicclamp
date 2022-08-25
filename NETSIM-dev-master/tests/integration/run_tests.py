#!/usr/bin/python3

import sys
from pathlib import Path
import subprocess

def main():
    #path to NETSIM dir
    if len(sys.argv) > 1:
        print("Using " + sys.argv[1] + " is the root NETSIM directory\n")
        root_path = Path(sys.argv[1]).absolute()
    else:
        print("Assuming the cwd is the root NETSIM directory\n")
        root_path = Path().cwd().absolute()

    inputs = ("test1.parameters",
              "test2.parameters",
              #"test3.parameters",
              "test4.parameters",
              "test5.parameters",
              "test6.parameters",
              "test7.parameters",
              "test8.parameters",
              "test9.parameters")    

    job = 0
    i = 0
    
    for test in inputs:
        print("Testing netsim with parameter file: " + test + " - Result: ", end="")
        r = subprocess.run([root_path / "netsim",
                            "-j " + str(job),
                            "-f", root_path / "tests" / test,
                            "-o", root_path / "data"
                            ],capture_output=True)
        
        if r.returncode != 0:
            print("Error in executing netsim \nstdout:", r.stdout.decode("utf-8"),
                  "\nstderr:", r.stderr.decode("utf-8"))
            exit()
            
        #get produced spike count
        for line in r.stdout.decode("utf-8").splitlines():
            if line.find("Total number of spikes:") != -1:
                produced_spikes = line.strip(" \n").split(" ")[4]
            
        
        #get true spike count    
        f = open(root_path / "tests" / test)
        for line in f:
            if line.find("Total Spikes:") != -1:
                expected_spikes = line.strip(" #\n").split(" ")[2]
        f.close()
        
        if produced_spikes == expected_spikes:
            print("\u001b[32m" + "Successful" + "\u001b[0m")
        else:
            print("\u001b[31m" + "Failed" + "\u001b[0m")

if __name__ == "__main__":
    main()
