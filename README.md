# Global-Sequence-Alignment
Purpose: Program determines the optimal global sequence alignment between two sequences. 

Directions on Compilation: 
    Compilation should follow the following format
    python3 align.py Match Mismatch Gap < input > output

    In this case, we are piping a FASTA-formatted file into stdin and sending
    all output to some output.txt file. Regarding the scoring scheme component,
    the user can pass in arguments (numbers) on the command line and they will
    be utilized to initialize the match, mismatch, and gap scores 
    respectively. However, if the user chooses not to do this, the program will
    still run, and it will use the default scores of match = 4, mismatch = -2,
    and gap = -2. 
