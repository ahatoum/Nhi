import subprocess

"""
Notes:

Don't forget hmmer python wrapper package (import hmmer) you downloaded. It has readymade methods and classes.

"""

def raxml():

    # subprocess.Popen(["ls", "-l"], cwd="/Users/baslan/raxml-ng_v1.0.2_macos_x86_64/")
    # list_files = subprocess.run(["ls", "-l"])
    # print("The exit code was: %d" % list_files.returncode)
    # list_files = subprocess.run(["raxml-ng", "-h"], cwd="/Users/baslan/raxml-ng_v1.0.2_macos_x86_64/")

    """

    Executing Shell Commands with Python
    https://stackabuse.com/executing-shell-commands-with-python/

    """

    # Works well
    subprocess.Popen(["./raxml-ng", "-h"], cwd="/Users/baslan/raxml-ng_v1.0.2_macos_x86_64/")

def hmmscan(path_input, path_output):

    # print(proc.communicate())
    # print("the commandline is {}".format(proc.args))
    # print(' '.join(command))

    """
    path_output = "works_yey.txt"  # '/Users/baslan/Biosoft/hmmer/works_yey.txt'
    path_input =  "sp2.fa"  #'/Users/baslan/Biosoft/hmmer/sp1.fa'

    path_output = "/Users/baslan/Biosoft/hmmer/works_yey.txt"
    path_input =  "/Users/baslan/PycharmProjects/NHI_proximity_To_Defences/neighbours/NZ_CP028295.1_WP_107133143.1.fa"
    """

    path_application = '/Users/baslan/Biosoft/hmmer/'

    command = ["./hmmscan",
               "--cut_ga",
               "--noali",
               "--acc",
               "--domtblout",
               path_output,
               "Pfam-A.hmm",
               path_input]

    proc = subprocess.Popen(command, cwd=path_application, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) # stdout & stderr arguments do the same job as '> /dev/null' if you were to type in command line.
    proc.wait()  # This makes sure python waits until the process is complete.


