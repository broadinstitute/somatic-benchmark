import sys
import os

def indelocatorLineMatches(line):
    "SOMATIC" in line and not "AUTOFILTER" in line

def strelkaLineMatches(line):
    not line.startswith("#")

def countIndels(path, lineMatches):
    try:
        indelCount = 0
        with open(path) as indelfile:
            for line in indelfile:
                if lineMatches(line):
                    indelCount+=1
        return indelCount
    except:
        print "Cannot read indels from {file}, setting count to -1.".format(file=path) 
        return -1

def findDirectories(baseDir, dirPrefix):
  return [baseDir+x for x in os.listdir(baseDir) if x.startswith(dirPrefix) and os.path.isdir("{base}/{x}".format(base=baseDir, x=x) ) ]

def outputToFile(results,outfile="indelocator_falsepositives_results.txt"):
    with open(outfile, 'w') as outfile:
        print "Writing to "+outfile.name
        outfile.write("Directory\tTool\tResult\n" )
        for dir,result in results.iteritems():
            outfile.write(formatOutputLine(dir,result))

def formatOutputLine(dir,result):
    results = ["{dir}\t{tool}\t{result}\n".format(dir=dir, tool=tool, result=result) for tool, result in result.iteritems()]
    return "".join(results)
      

def main(argv):
    if len(argv) <2:
        print "Usage <{thisScript}> <directory>".format(thisScript=argv[0]) 
        print "directory should contain all the generated indl_<lbl> subdirectories"
        quit(1)

    dir = argv[1]
    print "Looking in {dir}".format(dir=dir)

    indelocatorNames ={ "Raw":"sample.indels.txt",
                "N":"sample.n_filtered.indels.txt",
                "T":"sample.t_filtered.indels.txt",
                "Panel":"sample.filtered.panel_marked.indels.txt"}

    strelkaNames = {"Strelka": "final.indels.vcf"}
    res1 = gatherResults(argv[1], "Indelocator_", indelocatorNames, indelocatorLineMatches())
    res2 = gatherResults(argv[1], "Strelka_", strelkaNames, strelkaLineMatches())

    outputToFile(res1++res2)

def gatherResults(rootDir, prefix, filenames, lineMatcher):
    directories = findDirectories(rootDir, "Indelocator_")
    print "Found: " + ", ".join(directories)
    results = {}
    for directory in directories:
        dirName = directory.split("/")[-1]
        dirResults = {}
        for (lbl,file) in filenames.iteritems():
            path = "{dir}/{file}".format(dir=directory, file=file)
            indels = countIndels(path, indelocatorLineMatches)
            dirResults[lbl]=indels
        results[dirName] = dirResults
        print formatOutputLine(directory, dirResults)
    return results


if __name__ == "__main__":
       main(sys.argv)
