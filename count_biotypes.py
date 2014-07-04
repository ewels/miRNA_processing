#!/usr/bin/python
"""
small RNA pipeline
"""

from __future__ import print_function

import argparse
import datetime
import HTSeq
import numpy
import os
import re
import shlex
import shutil
import sys
import subprocess
import tempfile

from collections import defaultdict

# Import matplot lib but avoid default X environment
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt



def main(annotation_file, input_bam_list, biotype_flag, counts_file, feature_type, num_lines, quiet):
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
        (biotype_counts, biotype_lengths, aligned_reads, counts_string) = count_biotypes(fname, annotation_file, biotype_flag, feature_type, num_lines, quiet)
        if counts_file:
            # Save to file
            try:
                with open(counts_file, 'w') as fh:
                    print(counts_string, file=fh);
            except IOError as e:
                raise IOError(e)
        else:
            # Print to STDOUT
            print(counts_string, file=sys.stdout)
        
        # Plot graph
        plot_basename = os.path.splitext(os.path.basename(fname))[0]
        (bargraph_png, bargraph_pdf) = plot_bars(biotype_counts, aligned_reads, plot_basename, quiet)
            
    # Done!
    pass



def count_biotypes (aligned_bam, annotation_file, biotype_label='gene_type', count_feature_type='exon', number_lines=10000000, quiet=0):
    """
    Custom function that uses HTSeq core to analyse overlaps
    with annotation features of different biotypes.
    Returns a dict of biotypes with their counts, a dict
    of biotypes with lists of read lengths, the total number
    of aligned reads and a summary string ready for printing
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
    biotype_lengths = {}
    biotype_counts['no_overlap'] = 0
    biotype_counts['multiple_features'] = 0
    biotype_lengths['no_overlap'] = []
    biotype_lengths['multiple_features'] = []
    aligned_reads = 0
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
            if biotype_label not in feature.attr and biotype_label == 'gene_type':
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
                biotype_lengths[ feature.attr[biotype_label] ] = []
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
        print("\nReading BAM file (each dot is 1000000 lines, will stop at {}): ".format(number_lines), end='', file=sys.stderr)
    i = 0
    for alnmt in bamfile:
        i += 1
        if i > int(number_lines):
            i -= 1
            print("Reached {} lines in the aligned file, exiting..".format(number_lines), file=sys.stderr)
            break
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
            
            key = 'multiple_features'
            if len(iset) == 1:
                key = list(iset)[0]
            elif len(iset) == 0:
                key = 'no_overlap'
            biotype_counts[key] += 1
            biotype_lengths[key].append(alnmt.iv.length)
    
    if not quiet:
        print ("\n\n{} overlaps found from {} aligned reads ({} reads total)".format(aligned_reads-biotype_counts['no_overlap'], aligned_reads, i), file=sys.stderr)
        print ("{} reads had multiple feature overlaps".format(biotype_counts['multiple_features']), file=sys.stderr)
    
    
    # Make a string table out of the counts
    counts_string = 'Type\tRead Count\n'
    for biotype in sorted(biotype_counts, key=biotype_counts.get, reverse=True):
        if biotype_counts[biotype] == 0:
            continue
        counts_string += "{}\t{}{}".format(biotype, biotype_counts[biotype], os.linesep)
    
    # Return the counts
    return (biotype_counts, biotype_lengths, aligned_reads, counts_string)





def plot_bars (biotype_counts, total_reads, output_basename, quiet=0):
    """
    Plots bar graph of alignment biotypes using matplotlib pyplot
    Input: dict of biotype labels and associated counts
    Input: total number of reads (for percentage axis)
    Input: output fn
    Returns filenames of PNG and PDF graphs
    """
    
    # SET UP VARIABLES
    bar_width = 0.8
    plt_labels = []
    plt_values = []
    plt_percentages = []
    for biotype in sorted(biotype_counts, key=biotype_counts.get):
        if biotype_counts[biotype] == 0:
            continue
        plt_labels.append(biotype)
        plt_values.append(biotype_counts[biotype])
        plt_percentages.append((biotype_counts[biotype]/total_reads)*100)
    ypos = numpy.arange(len(plt_labels))
    
    # SET UP OBJECTS
    fig = plt.figure()
    axes = fig.add_subplot(111)
    
    # PLOT GRAPH
    plt.bar(1, bar_width, plt_values, ypos, align='center', orientation='horizontal', linewidth=0)
    
    # Y AXIS
    axes.set_yticks(ypos)
    axes.set_yticklabels(plt_labels)
    axes.tick_params(left=False, right=False)
    axes.set_ylim((1,len(plt_labels)))
    
    # X AXIS
    axes.set_xscale('log')
    axes.grid(True, zorder=0, which='both', axis='x', linestyle='-', color='#DEDEDE', linewidth=1)
    axes.set_axisbelow(True)
    
    # LABELS
    axes.set_xlabel('Number of Reads')
    axes.set_ylabel('Biotype')
    axes.set_title("Annotation Biotype Alignments\n{}".format(output_basename))
    axes.tick_params(axis='both', labelsize=8)
    # Give more room for the labels on the left and top
    plt.subplots_adjust(left=0.22)
    
    # SAVE OUTPUT
    png_fn = "{}.png".format(output_basename)
    pdf_fn = "{}.pdf".format(output_basename)
    print("Saving to {} and {}".format(png_fn, pdf_fn), file=sys.stderr)
    fig.savefig(png_fn)
    fig.savefig(pdf_fn)
    
    # Return the filenames
    return(png_fn, pdf_fn)


def plot_epic_histogram (biotype_lengths, output_basename, quiet=0):
    """
    Plot awesome histogram of read lengths, with bars broken up by feature
    biotype overlap
    """
    # http://matplotlib.org/1.2.1/examples/pylab_examples/histogram_demo_extended.html
    # http://matplotlib.org/examples/pylab_examples/bar_stacked.html
    # n, bins, patches = P.hist(x, 10, normed=1, histtype='bar', stacked=True)
    


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Count read overlaps with different biotypes.")
    # TODO should allow multiple reference genomes but then need to determine how annotation files and reference files are linked
    parser.add_argument("-g", "--genome-feature-file", dest="annotation_file", required=True,
                        help="GTF/GFF genome feature file to use for annotation (must match reference file)")
    parser.add_argument("-t", "--genome-feature", dest="feature_type", default='exon',
                        help="Type of annotation feature to count")
    parser.add_argument("-b", "--biotype-flat", dest="biotype_flag", default='gene_type',
                        help="GTF biotype flag (default = gene_type or *biotype*)")
    parser.add_argument("-c", "--counts-output-file", dest="counts_file",
                        help="Output filename for biotype counts (will print to STDOUT if not specified)")
    parser.add_argument("-n", "--num-lines", dest="num_lines", default=10000000,
                        help="Number of alignments to query")
    parser.add_argument("-q", "--quiet", dest="quiet",
                        help="Suppress status messages sent to STDERR")
    parser.add_argument("input_bam_list", metavar='<BAM file>', nargs="+",
                        help="List of input BAM filenames")

    kwargs = vars(parser.parse_args())
    
    main(**kwargs)
    