import glob, os, sys

for wig_fname in glob.glob(sys.argv[1]+'/*.wig'):
    cmnd = 'nice /home/fjacobs/miranda/wigToBigWig '+wig_fname' '+sys.argv[2]+wig_fname.replace('.wig', '.bigwig')
    print cmnd, os.system(cmnd)

