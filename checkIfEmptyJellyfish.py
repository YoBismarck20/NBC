
import os
import sys



if __name__ == "__main__":

    f = open("finalDebug","r")
    previousLine=''
    for line in f:
        theLine = line
        if theLine == "#" and previousLine == "#":
            print("ERROR")
            sys.exit()
        else:
            previousLine = theLine
    
