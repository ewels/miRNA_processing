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



def main(annotation_file, input_bam_list, biotype_flag, counts_file, feature_type, quiet):
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
        if not quiet:
            print("Processing {}".format(fname), file=sys.stderr)
        counts = count_biotypes(fname, annotation_file, biotype_flag, feature_type, quiet)
        if counts_file:
            # Save to file
            try:
                with open(counts_file, 'w') as fh:
                    print(counts, file=fh);
            except IOError as e:
                raise IOError(e)
        else:
            # Print to STDOUT
            print(counts, file=sys.stdout)
            
    # Done!
    pass



def count_biotypes (aligned_bam, annotation_file, biotype_label, count_feature_type='exon', quiet=0):
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
    # Help from http://www-huber.embl.de/users/anders/HTSeq/doc/tour.html#tour
    if not quiet:
        print("\nParsing annotation file {}".format(annotation_file), file=sys.stderr)
        print("Each dot is 100000 lines: ", end='', file=sys.stderr)
    selected_features = HTSeq.GenomicArrayOfSets( "auto", stranded=False )
    ignored_features = 0
    used_features = 0
    biotype_counts = {}
    aligned_reads = 0
    no_overlap_reads = 0
    multimap_reads = 0
    feature_type_counts = defaultdict(int)
    feature_type_biotype_counts = defaultdict(lambda: defaultdict(int))
    i = 0
    for feature in gtffile:
        i += 1
        if i % 100000 == 0 and not quiet:
            print(".", end='', file=sys.stderr)
        # Collect features and initialise biotype count objects
        if feature.type == count_feature_type:
            # See if we have another annotation that sounds like biotype
            # eg. Human ensembl calls it gene_biotype
            if biotype_label not in feature.attr and biotype_label == 'biotype':
                for attr in feature.attr.iterkeys():
                    if 'biotype' in attr:
                        if not quiet:
                            print("\nChanging biotype label from {} to {}".format(biotype_label, attr), file=sys.stderr)
                        biotype_label = attr
            
            # Initiate count object and add feature to selected_features set  
            if biotype_label in feature.attr:
                used_features += 1
                selected_features[ feature.iv ] += feature.attr[biotype_label]
                biotype_counts[ feature.attr[biotype_label] ] = 0
            else:
                ignored_features += 1
                
        # Collect general annotation stats
        feature_type_counts[feature.type] += 1
        if biotype_label in feature.attr:
            feature_type_biotype_counts[feature.type][feature.attr[biotype_label]] += 1
    
    if not quiet:
        print("\n\n{} features with biotype: {}".format(count_feature_type, used_features), file=sys.stderr)
        print("{} features without biotype: {}".format(count_feature_type, ignored_features), file=sys.stderr)
        print("{} biotypes to be counted: {}".format(count_feature_type, ', '.join(biotype_counts.keys())), file=sys.stderr)
        
        print("\nBiotype stats found for all feature types (using attribute '{}'):".format(biotype_label), file=sys.stderr)
        for ft in sorted(feature_type_biotype_counts.keys()):
            num_ft_bts = len(feature_type_biotype_counts[ft].keys())
            num_features = 0
            for c,d in feature_type_biotype_counts[ft].iteritems():
                num_features += d
            print("    {:20}\t{:4} biotypes\t{:6} labelled features".format(ft, num_ft_bts, num_features), file=sys.stderr)
    
    if(used_features == 0):
        raise Exception('No features have biotypes!')
    
    # Go through alignments, counting transcript biotypes
    if not quiet:
        print("\nReading BAM file (each dot is 1000000 lines): ", end='', file=sys.stderr)
    i = 0
    for alnmt in bamfile:
        i += 1
        if i % 1000000 == 0 and not quiet:
            print(".", end='', file=sys.stderr)
        if alnmt.aligned:
            aligned_reads += 1
            iset = None
            for iv2, step_set in selected_features[ alnmt.iv ].steps():
                if iset is None:
                    iset = step_set.copy()
                else:
                    iset.intersection_update( step_set )

            if len( iset ) == 1:
                biotype_counts[ list(iset)[0] ] += 1
                print("Name: >{}<".format(list(iset)[0]), file=sys.stderr)
            elif len(iset) == 0:
                no_overlap_reads += 1
            else:
                multimap_reads += 1
    
    if not quiet:
        print ("\n\n{} overlaps found from {} aligned reads ({} reads total)".format(aligned_reads-no_overlap_reads, aligned_reads, i), file=sys.stderr)
        print ("{} reads discarded due to multiple feature overlaps".format(multimap_reads), file=sys.stderr)
    
    
    # Print the counts
    if not quiet:
        print("\n\nType\tRead Count", file=sys.stderr)
    
    counts_string = ''
    for biotype in sorted(biotype_counts.keys()):
        counts_string += "{}\t{}{}".format(biotype, biotype_counts[biotype], os.linesep)
    
    # Return the counts
    return (counts_string)





if __name__ == "__main__":
    parser = argparse.ArgumentParser("Count read overlaps with different biotypes.")
    # TODO should allow multiple reference genomes but then need to determine how annotation files and reference files are linked
    parser.add_argument("-g", "--genome-feature-file", dest="annotation_file",
                        help="GTF/GFF genome feature file to use for annotation (must match reference file)")
    parser.add_argument("-f", "--genome-feature", dest="feature_type", default='exon',
                        help="Type of annotation feature to count")
    parser.add_argument("-b", "--biotype-flat", dest="biotype_flag", default='biotype',
                        help="GTF biotype flag (default = *biotype*)")
    parser.add_argument("-c", "--counts-output-file", dest="counts_file",
                        help="Output filename for biotype counts (will print to STDOUT if not specified)")
    parser.add_argument("-q", "--quiet", dest="quiet",
                        help="Suppress status messages sent to STDERR")
    parser.add_argument("input_bam_list", metavar='<BAM file>', nargs="+",
                        help="List of input BAM filenames")

    kwargs = vars(parser.parse_args())
    
    main(**kwargs)
    