#!/user/bin/env python3
""" Author: Lukas van den Heuvel Student number: 951218334040 Usage: 
n.a. Description: Subroutine template for command-line tools """
# Import statements
from sys import argv
import subprocess
import os.path

def run_program(infile1,infile2,infile3,output_filename):
    """Run the program ...
    output_filename -- the name the output should be saved to
    return: N/A, the function writes to a file
    """
    error = False
    # check if the input files exist, if not return an error exit()
    if os.path.exists(infile1):
        pass
    else:
        print("there is no such file: {}".format(infile1))
        error = True

    if os.path.exists(infile2):
        pass
    else:
        print("there is no such file: {}".format(infile2))
        error = True
    if os.path.exists(infile3):
        pass
    else:
        print("there is no such file: {}".format(infile3))
        error = True

    # check if the output file already exists, if not run the program
    if os.path.exists(output_filename):
        print("{0} already exists\n".format(output_filename))

        while question_loop == True:
            overwrite_question = input("Would you like to overwrite it? [y/n]")
            if overwrite_question == "y":
                check_call("rm {}".format(outfilename))
                question_loop = False
            elif overwrite_question == "n":
                error = True
            else:
                print("Please provide an answer as 'y' (yes) or 'n' (no)")
    if error == True:
        print("The pipeline was interupted due to an error in the ASSEMBLY")
        print("If you wish to restart your pipeline at the ASSEMBLY:")
        print("\tRun the SCRIPTNAME with ar argument '--step ASSEMBLY'\n") #ADD SCRIPT NAME!!!!!!!!!!
        exit()
    else:
        cmd_command = "cat {} {} {} > {}".\
        format(infile1,infile2,infile3,output_filename)
        exit_message = subprocess.check_call(cmd_command, shell=True)
        # feedback, not needed
        print("EXIT STATUS: {0}".format(exit_message))
        print("{0} was created\n".format(output_filename))
# main
if __name__ == '__main__':
    infile1 = "file1.txt"
    infile2 = "file2.txt"
    infile3 = "file3.txt"
    output_filename = "file.out"
    run_program(infile1,infile2,infile3,output_filename)
