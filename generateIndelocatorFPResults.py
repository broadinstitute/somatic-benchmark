import sys
import os

def countIndels(path):
    try:
        indelCount = 0
        with open(path) as indelfile:
            for line in indelfile:
                if "SOMATIC" in line and not "AUTOFILTER" in line:
                    indelCount+=1
        return indelCount
    except:
        print "Cannot read indels from {file}, setting count to -1.".format(file=path) 
        return -1

def findDirectories(baseDir):
  return [x for x in os.listdir(baseDir) if x.startswith("indl_") and os.path.isdir(x)] 

def outputToFile(results,outfile="indelocator_falsepositives_results.txt"):
    with open(outfile, 'w') as outfile:
        print "Writing to "+outfile.name
        outfile.write("Directory\tRaw\tN\tT\tPanel\n" )
        for dir,result in results.iteritems():
            outfile.write(formatOutputLine(dir,result))

def formatOutputLine(dir,result):
    return  "{dir}\t{raw}\t{n}\t{f}\t{panel}\n".format(dir=dir,raw=result["raw"],
                                                        n=result["N"],f=result['T'],
                                                        panel=result["Panel"])
      

def main(argv):
    if len(argv) <2:
        print "Usage <{thisScript}> <directory>".format(thisScript=argv[0]) 
        print "directory should contain all the generated indl_<lbl> subdirectories"
        quit(1)

    directories = findDirectories(argv[1]) 
    
    filenames ={ "raw":"indels.txt",
                "N":"n_filtered.indels.txt",
                "T":"t_filtered.indels.txt",
                "Panel":"filtered.panel_marked.indels.txt"}
    individual = "sample"
   
    results = {}
    for directory in directories:
        dirName = directory.split("/")[-1]
        dirResults = {}
        for (lbl,file) in filenames.iteritems():
            path = "{dir}/{indv}.{file}".format(dir=directory,indv=individual,file=file)
            indels = countIndels(path)
            dirResults[lbl]=indels
        results[dirName] = dirResults
        print formatOutputLine(directory, dirResults)
  
    outputToFile(results) 


if __name__ == "__main__":
       main(sys.argv)
