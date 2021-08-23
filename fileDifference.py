
import os
import sys


if __name__=="__main__":
    file1 = "filesRead"
    file2 = "fileOut"
    
    theSet = set()
    theSet2 = set()
    with open(file1,"r") as f:
        for line in f:
            theSet.add(line)
    with open(file2,"r") as f:
        for line in f:
            theArr = line.split("@")
            theSet2.add(theArr[1])

    print(len(theSet.difference(theSet2)))
    for stuff in theSet.difference(theSet2):
        print(stuff)
