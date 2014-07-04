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

import count_biotypes

from collections import defaultdict
from matplotlib import pyplot


def main(genomeref_file, annotation_file, id_attr, feature_type, mirbase_file, output_dir, num_cores, force_overwrite, keep_tmp, tmp_dir, quiet, input_fastq_list):
    """
    Add docstring here
    """
    
    # Sanity check - make sure input files exist
    for fname in input_fastq_list:
        if not os.path.isfile(fname):
            raise IOError("Fatal error - can't find input file {}".format(fname))
    
    # Load environment modules
    env_modules = ['bioinfo-tools','cutadapt','FastQC/0.11.2','bowtie','samtools', 'htseq']
    load_modules(env_modules)
    
    # Find the number of available cores
    num_cores = get_validated_cores(num_cores)
    
    # Set up directories
    run_directory = RunDirectory(output_dir, tmp_dir, keep_tmp, force_overwrite)
    
    # TODO Initialize logger -- see how Guillermo tends to do this -- need to write to local file, stderr, and logstash
    # logger = 
    
    # Merge and decompress input files
    working_fastq = merge_input_fastq_files(input_fastq_list, run_directory)
    
    # cutadapt
    trimmed_fastq = trim_fastq(working_fastq, run_directory)

    # fastqc
    fastqc_report = fastqc(trimmed_fastq, run_directory)

    # bowtie alignment
    aligned_bam = bowtie_align(trimmed_fastq, run_directory, genomeref_file, num_cores)
    
    # annotate with htseq-count
    if annotation_file and os.path.isfile(annotation_file):
        # Get top hits with HTSeq-counts on command line
        annotated_bam, htseq_counts_csv = htseq_counts(aligned_bam, run_directory, annotation_file, id_attr)
        # Make pie chart of biotype alignments
        (counts, piechart) = htseq_biotypes(aligned_bam, run_directory, annotation_file, feature_type)
    else:
        print("Error - annotation file not found, cannot calculate feature enrichment. " \
              "Use -g to specify a GTF/GFF file.\nSkipping feature enrichment step..\n\n",
               file=sys.stderr)
    
    # TODO - offer secondary annotation against miRBase co-ordinates?
    # http://www.mirbase.org/ftp.shtml
    
    # Align against miRBase
    miRBase_hairpin_aligned, miRBase_mature_aligned = miRBase_align(trimmed_fastq, run_directory, num_cores)
    
    # Sort and index the miRBase alignments
    hairpin_sorted, hairpin_index = bam_sort_index(miRBase_hairpin_aligned, run_directory)
    mature_sorted, mature_index = bam_sort_index(miRBase_mature_aligned, run_directory)
    
    # Produce counts of miRNA alignments
    hairpin_counts = samtools_count_alignments(hairpin_sorted, run_directory)
    mature_counts = samtools_count_alignments(mature_sorted, run_directory)
    
    # visualizations

    # remove or save tmp directory
    run_directory.remove_tmp_dir
    
    pass

def load_modules(modules):
    """
    Takes a list of environment modules to load (in order) and
    loads them using modulecmd python load
    Returns True
    """
    # Module loading is normally controlled by a bash function
    # As well as the modulecmd bash which is used in .bashrc, there's also
    # a modulecmd python which allows us to use modules from within python
    # UPPMAX support staff didn't seem to know this existed, so use with caution
    for mod in modules:
        p = subprocess.Popen("modulecmd python load "+mod,
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        stdout,stderr = p.communicate()
        try:
            exec stdout
        except:
            print("Error loading modules:", stderr, file=sys.stderr)
    return True
          
        
def get_validated_cores(num_cores):
    sys_cores =     os.getenv('SLURM_CPUS_ON_NODE') \
                or  subprocess.check_output(['nproc', '--all']) \
                or  1
    sys_cores = int(sys_cores)

    if not num_cores:
        num_cores = sys_cores
    if num_cores > sys_cores:
        print(  "Requested number of cores ({num_cores}) greater than number of system cores ({sys_cores}); " \
                "using {sys_cores} instead.".format(**locals()), file=sys.stderr)
        num_cores = sys_cores
    if num_cores < 1:
        print(  "Requested number of cores ({num_cores}) must be a postive integer; " \
                "using 1 instead.".format(**locals()), file=sys.stderr)
        num_cores = 1
    return num_cores

def is_compressed(input_file):
    """
    Checks to see if the file is gzip compressed, returns T/F
    """
    # TODO add bzip2 support

    # Challenge accepted.. (file seems to be much faster than gzip -d -t)
    cmd = shlex.split("file {}".format(input_file))
    try:
        file_output = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        out, err = file_output.communicate()
        if 'ASCII text' in out:
            return False
        elif 'compressed' in out:
            return True
        else:
            print ("Didn't recognise type of input file: {}".format(out), file=sys.stderr)
            raise
    except subprocess.CalledProcessError, e:
        print ("Couldn't check type of input file: {}".format(e.output), file=sys.stderr)
        raise

def decompress_file(input_file, output_dir=os.getcwd(), return_pipe=False):
    """
    Decompresses an input file in gzip format.
    Returns the abspath to the decompressed output file
    or a PIPE if requested.
    """
    # TODO add bzip2 support

    if is_compressed(input_file):
        if return_pipe:
            cmd = shlex.split("gzip -d -c {0}".format(input_file))
            return subprocess.Popen(cmd, stdout=subprocess.PIPE)
        else:
            # Remove .gz, .gzip, .bz2 extensions if present
            basename = strip_gzip_ext(os.path.basename(input_file))
            output_file = os.path.join(output_dir, basename)
            print ("Unzipping {} to {}..".format(input_file, output_file), file=sys.stderr)
            cmd = shlex.split("gzip -d -c {0}".format(input_file))
            with open(output_file, 'w') as f:
                subprocess.Popen(cmd, stdout=f)
            return output_file
    else:
        if return_pipe:
            cmd = shlex.split("cat {}".format(input_file))
            return subprocess.Popen(cmd, stdout=subprocess.PIPE)
        else:
            return os.path.realpath(input_file)

def strip_gzip_ext(file_name):
    """
    Returns the file without the .gz or .gzip suffix, if present
    """
    # TODO add bzip2 support

    base, ext = os.path.splitext(file_name)
    # if re.match('\.gz(ip)?$|.bz2', ".bz2", ext):
    if re.match('\.gz(ip)?$|.bz2', ext):
        file_name = base
    return file_name


def merge_input_fastq_files(input_fastq_list, run_directory):
    """
    Merge multiple fastq files into one fastq file.
    Returns the path to the final merged fastq file,
    or the original file in the case of just one input file.
    """
    if len(input_fastq_list) is 1:
        fastq_file = input_fastq_list[0]
        # decompress if compressed
        return decompress_file(fastq_file, output_dir=run_directory.tmp_dir)
    else:
        merged_filename = "MERGED_{}".format(strip_gzip_ext(os.path.basename(input_fastq_list[0])))
        merged_filepath = os.path.join(run_directory.tmp_dir, merged_filename)
        with open(merged_filepath, 'w') as output_file:
            for fastq_file in input_fastq_list:
                stream = decompress_file(fastq_file, return_pipe=True)
                output_file.write(stream.stdout.read())
        return merged_filepath

def trim_fastq(fq_input, run_directory, min_qual=10, min_match=3, min_length=18, adapter="TGGAATTCTCGGGTGCCAAGG"):
    """
    Run Cutadapt with a FastQ input file
    """
    fq_output_fn = "{}_trimmed.fq".format(os.path.splitext(os.path.basename(fq_input))[0])
    fq_output = os.path.join(run_directory.tmp_dir, fq_output_fn)
    
    # Put the command together
    cmd = shlex.split("cutadapt -f fastq -a {adapter} -q {min_qual} " \
                      "--match-read-wildcards -O {min_match} -m {min_length} " \
                      "-o {fq_output} {fq_input}".format(**locals()))
    # Status message
    timestamp = datetime.date.strftime(datetime.datetime.now(), format="%H:%M:%S %d/%m/%Y")
    print("\nRunning cutadapt. Started at {0}\nCommand: {1}\n".format(timestamp,
                                                    " ".join(cmd)), file=sys.stderr)
    
    # Run the command
    try:
        subprocess.check_call(cmd)
        return fq_output
    except subprocess.CalledProcessError:
        return exit(1)



def fastqc(fq_input, run_directory):
    """
    Run FastQC with a FastQ input file.
    Returns the directory path containing the FastQC output files
    """
    outdir = os.path.join(run_directory.output_dir, "FastQC")
    
    # Make this directory if it doesn't already exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    
    # Put the command together
    # When FastQC is updated on UPPMAX, can use -q to prevent printing out %ages
    cmd = shlex.split("fastqc -o {outdir} {fq_input}".format(**locals()))
    # Status message
    timestamp = datetime.date.strftime(datetime.datetime.now(), format="%H:%M:%S %d/%m/%Y")
    print("\nRunning FastQC. Started at {0}\nCommand: {1}\n".format(timestamp,
                                                    " ".join(cmd)), file=sys.stderr)
    
    # Run the command
    try:
        subprocess.check_call(cmd)
        return outdir
    except subprocess.CalledProcessError:
        return False


def bowtie_align(fq_input, run_directory, genomeref_file, num_cores):
    """
    Run Bowtie 1 with a FastQ input file and create aligned BAM file.
    Uses parameters best suited to miRNA alignment.
    Returns the path to the aligned BAM file.
    """
    fq_input = os.path.realpath(fq_input)
    genomeref_file = os.path.realpath(genomeref_file)
    bam_output_fn = "{}_aligned.bam".format(os.path.splitext(os.path.basename(fq_input))[0])
    bam_output = os.path.join(run_directory.output_dir, bam_output_fn)
    
    # Put the bowtie 1 command together...
    # -p {num_cores}        Number of cores to use for alignment
    # -t                    Print out timestamps
    # -n 0                  Do not allow any mismatches in seed
    # -l 15                 Seed length of 15bp
    # -e 99999              Disable max sum of mismatch quals across alignment
    # -k 200                Report up to 200 good alignments per read
    # --best                Report hits from best stratum
    # -S                    SAM output
    # --chunkmbs 2048       Up the max mb RAM from 64mb to 2gigs    
    b_cmd = shlex.split("bowtie -p {num_cores} -t -n 0 -l 15 -e 99999 -k 200 " \
                        "--best -S --chunkmbs 2048 {genomeref_file} {fq_input}" \
                        .format(**locals()))
    
    # Write the samtools command to convert the SAM to a BAM
    # We'll write the stdout from this command to a file using Popen
    samtools_cmd = shlex.split("samtools view -bSh -")
    
    # Status message
    timestamp = datetime.date.strftime(datetime.datetime.now(), format="%H:%M:%S %d/%m/%Y")
    print("\nRunning Bowtie alignment. Started at {0}\nBowtie Command: {1}\n" \
          "Samtools Command: {2}\n".format(timestamp, " ".join(b_cmd),
                                      " ".join(samtools_cmd)), file=sys.stderr)
    
    # Run the command
    try:
        with open(bam_output, 'w') as fh:
            bowtie = subprocess.Popen(b_cmd, stdout=subprocess.PIPE)
            samtools = subprocess.Popen(samtools_cmd, stdin=bowtie.stdout, stdout=fh).communicate()
            bowtie.stdout.close()
        return bam_output
        
    except subprocess.CalledProcessError:
        return exit(1)
        
        
    
def htseq_counts(aligned_bam, run_directory, annotation_file, id_attr='gene_id'):
    """
    Run htseq-count on the command line to create an exon counts file sorted by count
    Returns paths to the annotated BAM file and counts file, as a list
    """
    # Set up filenames
    aligned_bam = os.path.realpath(aligned_bam)
    annotation_file = os.path.realpath(annotation_file)
    annotated_bam_fn = "{}_annotated.bam".format(os.path.splitext(os.path.basename(aligned_bam))[0])
    annotated_bam = os.path.join(run_directory.output_dir, annotated_bam_fn)
    htseq_counts_fn = "{}_counts.txt".format(os.path.splitext(os.path.basename(aligned_bam))[0])
    htseq_counts = os.path.join(run_directory.output_dir, htseq_counts_fn)
    
    # Write the commands
    samtools_view_cmd = shlex.split("samtools view -h {}".format(aligned_bam))
    htseq_cmd = shlex.split("htseq-count -o {} -t exon -s no -q -i '{}' - {}" \
                        .format(annotated_bam, annotation_file, id_attr))
    sort_cmd = shlex.split("sort -n -k 2 -r")
    
    # Pipe, baby, pipe
    counts_output = ""
    try:
        p1 = subprocess.Popen(samtools_view_cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(htseq_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)
        p1.stdout.close()
        p3 = subprocess.Popen(sort_cmd, stdin=p2.stdout, stdout=subprocess.PIPE)
        p2.stdout.close()
        counts_output = p3.stdout.read()
    except subprocess.CalledProcessError:
        return exit(1)
    
    # Write to file, stripping features with zero counts (unless they're __special)
    try:
        with open(htseq_counts, 'w') as fh:
            print("Name\tRead Count\n", file=fh);
            for counts_line in counts_output.split(os.linesep):
                try:
                    [name, readcount] = counts_line.split('\t')
                except IndexError, e:
                    raise IndexError(e.args)
                if (readcount > 0 or name[:2] == '__'):
                    print(counts_line+os.linesep, file=fh);
    except IOError as e:
        raise IOError(e)
    
    # Return with filenames
    return (annotated_bam, htseq_counts)


def htseq_biotypes (bam_input, run_directory, annotation_file, feature_type='exon'):
    """
    Use count_biotypes.py to count read overlaps with different biotypes.
    """
    # Set up filenames
    bam_input = os.path.realpath(bam_input)
    annotation_file = os.path.realpath(annotation_file)
    counts_fn = "{}_counts.txt".format(os.path.splitext(os.path.basename(fname))[0])
    counts_path = os.path.join(run_directory.output_dir, counts_fn)
    
    # Get the counts
    counts = count_biotypes.count_biotypes(bam_input, annotation_file, count_feature_type, feature_type, quiet)
    
    # Save counts to file
    try:
        with open(counts_path, 'w') as fh:
            print(counts, file=fh);
    except IOError as e:
        raise IOError(e)
    
    
def miRBase_align(fq_input, run_directory, num_cores):
    """
    Run Bowtie 1 with a FastQ input file and creates an aligned BAM file.
    Aligns against miRBase hairpin and mature references
    Uses parameters best suited to miRNA alignment.
    Returns the path to the two aligned BAM files (hairpin, mature).
    """
    fq_input = os.path.realpath(fq_input)
    miRBase_dir = os.path.dirname(os.path.realpath(__file__)) + '/miRBase/'
    miRBase_hairpin = miRBase_dir + 'hairpin'
    miRBase_mature = miRBase_dir + 'mature'
    hairpin_output_fn = "{}_hairpin_miRNA_aligned.bam".format(os.path.splitext(os.path.basename(fq_input))[0])
    hairpin_output = os.path.join(run_directory.tmp_dirs, hairpin_output_fn)
    mature_output_fn = "{}_mature_miRNA_aligned.bam".format(os.path.splitext(os.path.basename(fq_input))[0])
    mature_output = os.path.join(run_directory.tmp_dirs, mature_output_fn)
    
    # Alignment commands
    h_cmd = shlex.split("bowtie -p {num_cores} -t -n 0 -l 15 -e 99999 -k 200 " \
                        "--best -S --chunkmbs 2048 {miRBase_hairpin} {fq_input}" \
                        .format(**locals()))
    
    m_cmd = shlex.split("bowtie -p {num_cores} -t -n 0 -l 15 -e 99999 -k 200 " \
                        "--best -S --chunkmbs 2048 {miRBase_mature} {fq_input}" \
                        .format(**locals()))
    
    # Samtools command - **ignore unmapped reads**
    samtools_cmd = shlex.split("samtools view -bS -F 0x4 -")
    
    # Status message
    timestamp = datetime.date.strftime(datetime.datetime.now(), format="%H:%M:%S %d/%m/%Y")
    print("\nRunning miRBase alignments. Started at {0}\nHairpin Alignment Command: {1}\n" \
          "Mature Alignment Command: {2}\nSamtools Command: {3}\n".format(timestamp, " ".join(h_cmd),
                                      " ".join(m_cmd), " ".join(samtools_cmd)), file=sys.stderr)
    
    # Run the hairpin alignment
    try:
        with open(hairpin_output, 'w') as fh:
            bowtie = subprocess.Popen(h_cmd, stdout=subprocess.PIPE)
            samtools = subprocess.Popen(samtools_cmd, stdin=bowtie.stdout, stdout=fh).communicate()
            bowtie.stdout.close()
    except subprocess.CalledProcessError:
        return exit(1)
    
    # Run the mature alignment
    try:
        with open(mature_output, 'w') as fh:
            bowtie = subprocess.Popen(m_cmd, stdout=subprocess.PIPE)
            samtools = subprocess.Popen(samtools_cmd, stdin=bowtie.stdout, stdout=fh).communicate()
            bowtie.stdout.close()
    except subprocess.CalledProcessError:
        return exit(1)
    
    # Return with filenames
    return (hairpin_output, mature_output)


def bam_sort_index (input_bam, run_directory):
    """
    Run samtools sort and then samtools index on a BAM file.
    Returns two filenames - sorted BAM and index file.
    """
    # Set up filenames
    input_bam = os.path.realpath(input_bam)
    input_bam_basename = os.path.splitext(os.path.basename(input_bam))[0]
    sorted_bam_fn = "{}_sorted".format(input_bam_basename)
    bam_index_fn = "{}_sorted.bai".format(input_bam_basename)
    sorted_bam = os.path.join(run_directory.output_dir, sorted_bam_fn)
    bam_index = os.path.join(run_directory.output_dir, bam_index_fn)
    
    # Write sorting command
    samtools_sort_cmd = shlex.split("samtools sort {0} {1}" \
                        .format(input_bam, sorted_bam))
    
    # Write index command
    samtools_index_cmd = shlex.split("samtools index {0}.bam {1}" \
                        .format(sorted_bam, bam_index))
    
    # Status message
    # print("\nSorting and indexing BAM file\n", file=sys.stderr)
    # print("Sorting Command: {}\n".format(" ".join(samtools_sort_cmd), file=sys.stderr)
    # print("Indexing Command: {}\n".format(" ".join(samtools_index_cmd), file=sys.stderr)
                                      
    # Run the indexing
    try:
        subprocess.check_call(samtools_sort_cmd)
        # Run the sorting
        try:
            subprocess.check_call(samtools_index_cmd)
            return ("{}.bam".format(sorted_bam), bam_index)
        except subprocess.CalledProcessError:
            return False
    except subprocess.CalledProcessError:
        return False
        
        
def samtools_count_alignments (input_bam, run_directory):
    """
    Run samtools idxstats to count alignments in a BAM file.
    Requires a *sorted* and *indexed* BAM file.
    Writes an output file and returns filename
    """
    # Set up filenames
    input_bam = os.path.realpath(input_bam)
    input_bam_basename = os.path.splitext(os.path.basename(input_bam))[0]
    counts_fn = "{}_counts.txt".format(input_bam_basename)
    counts_path = os.path.join(run_directory.output_dir, counts_fn)
    
    # Write commands
    samtools_cmd = shlex.split("samtools idxstats {}".format(input_bam))
    sort_cmd = shlex.split("sort -r -k3,3")
    
    # Output
    counts_output = ""
    
    # Run the commands, save output in a string
    try:
        samtools_p = subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE)
        sort_p = subprocess.Popen(sort_cmd, stdin=samtools_p.stdout, stdout=subprocess.PIPE)
        samtools_p.stdout.close()
        counts_output = sort_p.communicate()[0]
    except subprocess.CalledProcessError:
        return False
    
    # Print a header and strip any lines with a count of 0
    try:
        with open(counts_path, 'w') as fh:
            print("Name\tLength\tNum Mapped\tNum Unmapped\n", file=fh);
            for counts_line in counts_output.split(os.linesep):
                try:
                    readcount = counts_line.split('\t')[2]
                except IndexError, e:
                    raise IndexError(e.args)
                if (readcount > 0):
                    print(counts_line+os.linesep, file=fh);
            
    except IOError, e:
        raise IOError(e.args)
    
    # Return counts filename
    return counts_path


class RunDirectory(object):
    """
    Keeps track of the various directories used for work and output.
    """
    # TODO I wonder if this is worth making into a context manager
    #       e.g. create tmp on entry, remove on exit

    def __init__(self, output_dir, tmp_dir, keep_tmp, force_overwrite):
        self.output_dir         = self.create_output_dir(output_dir)
        self.tmp_dir            = self.create_tmp_dir(tmp_dir)
        self.force_overwrite    = force_overwrite
        self.keep_tmp           = keep_tmp
        self.res_dirs           = []
        self.tmp_dirs           = []

    def create_output_dir(self, output_dir):
        """
        Create the output directory passed in by the user
        or one following the format ./"smRNA_run_(datetime)/"
        Returns the absolute path to the output directory.
        """
        if not output_dir:
            output_dir = os.path.join(  os.getcwd(),
                                        "smRNA_run_{}".format(datetime.date.strftime(
                                        datetime.datetime.now(), format="%Y%m%d_%H-%M-%S")))
        else:
            output_dir = os.path.realpath(output_dir)

        try:
            print('Creating output directory "{}"'.format(output_dir), file=sys.stderr)
            os.makedirs(output_dir)
        except OSError as e:
            if e.errno is 17:
            # Output directory already exists
                pass
        self.output_dir = output_dir
        return self.output_dir


    def create_tmp_dir(self, tmp_dir):
        """
        Create a temporary working directory; use the value passed by the user
        else system tmp if determinable else output directory.
        Returns the absolute path to the tmp directory.
        """
        # if keep_tmp, still write to tmp_dir so as not to hammer network drives, etc. (?)
        if not tmp_dir:
            # try using environment vars to locate system tmp
            tmp_dir = os.getenv('TMPDIR') or os.getenv('SNIC_TMP') or self.output_dir
        
        # Use the tempfile package to create a temporary directory
        tmp_dir = tempfile.mkdtemp(prefix='miRNAtmp_', dir=tmp_dir)
        print('Creating tmp directory "{}"'.format(tmp_dir), file=sys.stderr)

        self.tmp_dir = tmp_dir
        return self.tmp_dir
        
    # PHIL NOTE - what does this do?
    def create_dir(self, dir_name, in_tmp=False):
        """
        Create a new directory within the working or tmp directory.
        """
        if in_tmp:
            full_dir_path = os.path.join(self.tmp_dir, dir_name)
        else:
            full_dir_path = os.path.join(self.output_dir, dir_name)

        print('Creating directory "{}"'.format(dir_name), file=sys.stderr)
        os.makedirs(full_dir_path)

        if in_tmp:
            self.tmp_dirs.append(full_dir_path)
        else:
            self.dirs.append(full_dir_path)

        return full_dir_path


    def remove_tmp_dir(self):
        """
        Remove the tmp directory and its contents, unless -k / --keep-tmp
        is specified, in which case move the tmp directory into the
        working directory. Returns T/F
        """
        if self.keep_tmp:
            output_dir_tmp = os.path.join(self.output_dir, 'tmp')
            print('--keep-tmp specified. Moving temporary directory {0} to {1}' \
                    .format(self.tmp_dir, output_dir_tmp), file=sys.stderr)
            try:
                shutil.move(self.tmp_dir, output_dir_tmp)
                return True
            except:
                return False
        else:
            print('Removing temporary directory {}'.format(self.tmp_dir), file=sys.stderr)
            try:
                shutil.rmtree(self.tmp_dir)
                return True
            except:
                return False
            
        
        
          

if __name__ == "__main__":
    parser = argparse.ArgumentParser("Execute the small RNA pipeline.")
    # TODO should allow multiple reference genomes but then need to determine how annotation files and reference files are linked
    parser.add_argument("-r", "--genome-reference-file", dest="genomeref_file", required=True,
                        help="The genome reference file against which to align.")
    parser.add_argument("-g", "--genome-feature-file", dest="annotation_file",
                        help="GTF/GFF genome feature file to use for annotation (must match reference file).")
    parser.add_argument("-i", "--genome-id-attr", dest="id_attr", default='gene_id',
                        help="Annotation ID to use for HTSeq Counts")
    parser.add_argument("-y", "--genome-feature", dest="feature_type", default='exon',
                        help="Type of annotation feature to count with HTSeq")
    parser.add_argument("-m", "--mirbase-file",
                        help="The miRBase reference file.")
    parser.add_argument("-o", "--output-dir",
                        help="The output directory.")
    parser.add_argument("-n", "--num-cores", type=int,
                        help="The number of cores to use for alignment.")
    parser.add_argument("-t", "--tmp-dir",
                        help="Optionally specify a temporary directory.")
    parser.add_argument("-f", "--force-overwrite", action="store_true", default=False,
                        help="Force overwrite of existing files.")
    parser.add_argument("-k", "--keep-tmp", action="store_true", default=False,
                        help="Keep temporary files after processing.")
    parser.add_argument("-q", "--quiet", dest="quiet", default='0',
                        help="Suppress status messages sent to STDERR")
    parser.add_argument("input_fastq_list",  metavar='<BAM file>', nargs="+",
                        help="List of input FastQ filenames (can be gzipped)")

    kwargs = vars(parser.parse_args())
    
    main(**kwargs)
