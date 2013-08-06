import sys
import os

def countIndels(path):
    try:
        indelCount = 0
        with open(path) as indelfile:
            for line in indelfile:
                if not line.startswith("#"):
                    indelCount+=1
        return indelCount
    except:
        print "Cannot read indels from {file}, setting count to -1.".format(file=path) 
        return -1

def findDirectories(baseDir):
  return [baseDir+x for x in os.listdir(baseDir) if x.startswith("Strelka") and os.path.isdir("{base}/{x}".format(base=baseDir, x=x) ) ] 

def outputToFile(results,outfile="strelka_falsepositives_results.txt"):
    with open(outfile, 'w') as outfile:
        print "Writing to "+outfile.name
        outfile.write("Directory\tStrelka\n" )
        for dir,result in results.iteritems():
            outfile.write(formatOutputLine(dir,result))

def formatOutputLine(dir,result):
    return  "{dir}\t{raw}\n".format(dir=dir,raw=result["raw"])
      

def main(argv):
    if len(argv) <2:
        print "Usage <{thisScript}> <directory>".format(thisScript=argv[0]) 
        print "directory should contain all the generated strelka_<lbl> subdirectories"
        quit(1)

    dir = argv[1]
    print "Looking in {dir}".format(dir=dir)
    directories = findDirectories(argv[1]) 
    print "Found: " + ", ".join(directories)
     
    filenames ={ "raw":"final.indels.vcf"}
   
    results = {}
    for directory in directories:
        dirName = directory.split("/")[-1]
        dirResults = {}
        for (lbl,file) in filenames.iteritems():
            path = "{dir}/{file}".format(dir=directory,file=file)
            indels = countIndels(path)
            dirResults[lbl]=indels
        results[dirName] = dirResults
        print formatOutputLine(directory, dirResults)
  
    outputToFile(results) 


if __name__ == "__main__":
       main(sys.argv)
