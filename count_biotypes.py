#!/usr/bin/python
"""
small RNA pipeline
"""

from __future__ import print_function

import argparse
import datetime
import HTSeq
import os
import re
import shlex
import shutil
import sys
import subprocess
import tempfile

from collections import defaultdict
from matplotlib import pyplot



def main(annotation_file, input_bam_list):
    """
    Count the biotypes
    """
    
    # Sanity check - make sure input files exist
    if annotation_file:
        if not os.path.isfile(annotation_file):
            raise IOError("Fatal error - can't find annotation file {}".format(annotation_file))
    else:
        raise IOError("Fatal error - annotation file not specified")
    for fname in input_bam_list:
        if not os.path.isfile(fname):
            raise IOError("Fatal error - can't find input file {}".format(fname))
    
    # Process files
    for fname in input_bam_list:
        counts = count_biotypes(fname, annotation_file)
        # Save to file
        counts_fn = "{}_counts.txt".format(os.path.splitext(os.path.basename(fname))[0])
        try:
            with open(counts_fn, 'w') as fh:
                print(counts_fn, file=fh);
        except IOError as e:
            raise IOError(e)
    
    # Done!
    pass



def count_biotypes (aligned_bam, annotation_file):
    """
    Custom function that uses HTSeq core to calculate percentage alignments
    across different annotation gene biotypes.
    Returns filename to pie chart graphic.
    """
    # Set up filenames
    aligned_bam = os.path.realpath(aligned_bam)
    annotation_file = os.path.realpath(annotation_file)
    
    # Create file objects
    bamfile = HTSeq.BAM_Reader( aligned_bam )
    gtffile = HTSeq.GFF_Reader( annotation_file )
            
    # Go through annotation
    biotype_feature_type = 'transcript'
    selected_features = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    ignored_features = 0
    biotype_counts = {}
    feature_types = {}
    for feature in gtffile:
        # Collect features and initialise biotype count objects
        if feature.type == biotype_feature_type:
            if feature.attr['biotype'] is None:
                ignored_features += 1
            else:
                selected_features[ feature.iv ] += feature.attr['biotype']
                biotype_counts[ feature.attr['biotype'] ] = 0
                
        # Collect general annotation stats
        feature_types[feature.type]['count'] += 1
        if feature.attr['biotype'] is not None:
            feature_types[feature.type]['biotype_count'] += 1
            feature_types[feature.type]['biotypes'][feature.attr['biotype']] += 1
    
    # Go through alignments, counting transcript biotypes
    for alnmt in bamfile:
        if alnmt.aligned:
            iset = None
            for iv2, step_set in selected_features[ alnmt.iv ].steps():
                if iset is None:
                    iset = step_set.copy()
                else:
                    iset.intersection_update( step_set )
            if len( iset ) == 1:
                biotype_counts[ list(iset)[0] ] += 1
    
    # Print the counts
    counts_string = ''
    for biotype in sorted( biotype_counts.keys() ):
        counts_string += biotype, counts[biotype]
    
    # Return the counts
    return (counts_string)





if __name__ == "__main__":
    parser = argparse.ArgumentParser("Execute the small RNA pipeline.")
    # TODO should allow multiple reference genomes but then need to determine how annotation files and reference files are linked
    parser.add_argument("-g", "--genome-feature-file", dest="annotation_file",
                        help="GTF/GFF genome feature file to use for annotation (must match reference file).")
    parser.add_argument("input_bam_list", nargs="+")

    kwargs = vars(parser.parse_args())
    
    main(**kwargs)
    