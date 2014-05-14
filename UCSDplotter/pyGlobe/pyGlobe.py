import analysis
from optparse import OptionParser
import cProfile

def main(options):
    ana = analysis.Analysis(options)
    ana.Run()

if __name__ == "__main__":  
    parser = OptionParser()
    parser.add_option("-b", "--blind", default=False, action="store_true", help="Do not plot data")
    parser.add_option("-i", "--input", default="~/sslept_v3.root", help="Input filename")
    parser.add_option("-o", "--output", default="output_sslept.root", help="Output filename")
    parser.add_option("-I", "--inputfile", default="inputfiles.dat", help="Original inputfiles.dat")
    parser.add_option("-P", "--plotvariables", default="plotvariables.dat", help="Plot definition file")
    parser.add_option("-t", "--treename", default="opttree", help="Name of the tree to process")
    parser.add_option("-m", "--mode", default="highestSumpt", help="Selection mode (highestSumpt, allPairs...)")
    parser.add_option("-d", "--debug", action="store_true", default=False, help="Activate profiling of the program.")
    (options, arg) = parser.parse_args()
        
    if (options.debug):
        cProfile.run("main(options)")
    else:
        main(options)
