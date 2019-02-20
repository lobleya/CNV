#! /usr/bin/env python3

import pysam as ps
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-q", "--qual", dest="qual",
                  help="filter by read quality",
                  default=30,
                  metavar="float")

parser.add_option("-b", "--bam",
                  dest="bam",
                  nargs=2,
                  help="bam file ",
                  default="test.bam",metavar="file")

parser.add_option("-f", "--frag",
                  dest="frag",
                  help="fragment size filter",
                  default=150,metavar="float")

(options, args) = parser.parse_args()




samfile = ps.AlignmentFile(bam, "rb")
sorted_sam=samfile + ".sorted"


ps.sort("-n -o ",
         sorted_sam,
        " -T . ",
         samfile)
filtered= samfile + "filtered"
samfile = ps.AlignmentFile(sorted_sam, "rb")
paired  = ps.AlignmentFile(filtered, "wb",
                           template=samfile)


for read in ps.fetch() :
    if( read.is_proper_pair and abs(read.tlen) > frag and read.flag("0x256") and read.mapping_quality > qual ):
        paired.write(read)

paired.close()
paired.build()
samfile.close()


exit

#--- only print reads with flag 256 and with frag size > frag

#PE reads                R1--------->                    <---------R2

#--- exit planet dust !!!
