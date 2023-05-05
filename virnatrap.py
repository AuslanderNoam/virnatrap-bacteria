"""
Model classes and functions to identify viral reads and assemble viral contigs from input fastq
"""
# Imports --------------------------------------------------------------------------------------------------------------
import os
import re
import sys
from multiprocessing import Pool,freeze_support
import numpy as np
import random
from ctypes import *
from tensorflow import get_logger, autograph
from tensorflow.keras.models import load_model
#from tensorflow.keras.preprocessing.sequence import pad_sequences
import glob
import argparse
#import yep
from datetime import datetime


DESCRIPTION = "Extract viral contigs from a directory with unmapped RNAseq reads fastq files and saves a file with contigs for each fastq in an output directory"

os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'


get_logger().setLevel('ERROR')
autograph.set_verbosity(0)


flatten = lambda l: [item for sublist in l for item in sublist]

PWD = os.getcwd()

# Constants ------------------------------------------------------------------------------------------------------------
DEFAULT_NUC_ORDER = {y: x for x, y in enumerate(["A", "C", "G", "T"])}
NUCLEOTIDES = [x for x in DEFAULT_NUC_ORDER.keys()]
SEGMENT_LENGTH = 76
SEARCHSUBLEN = 24

# Functions ------------------------------------------------------------------------------------------------------------

def random_base(foo=""):
    """
    Generate a random base.
    :return: Random base.
    """
    return random.choice(NUCLEOTIDES)


def handle_N(sequence):
    """
    Handle non ATGCs.
    :param sequence: String input.
    :return: String output (only ATCGs), with randomly assigned bp to non-ATGCs.
    """
    ret = re.sub('N', random_base, sequence)
    return ret



def encode_sequence(sequence, nuc_order=None):
    """
    Encode a sequence to integers for use in Python LSTM model.
    :param sequence: Sequence to encode.
    :param nuc_order: Order of nucleotides for encoding.
    :return: Encoded sequence as integers.
    """
    if nuc_order is None:
        nuc_order = DEFAULT_NUC_ORDER

    #sequence = sequence.upper()[:SEGMENT_LENGTH]
    #accepted_nucleotides = "ACGT"

    #assert re.match('^[{}]+$'.format(accepted_nucleotides), sequence) is not None, \
    #    "Only {} allowed".format(accepted_nucleotides)

    encoded_seq = [nuc_order[x] for x in sequence]
    #encoded_seq = np.array([encoded_seq])

    return encoded_seq


def encode_sequences(sequences, nuc_order=None, segment_length=SEGMENT_LENGTH):
    """
    Encode a sequence to integers for use in model.
    :param sequences: List of sequences to encode.
    :param nuc_order: Order of nucleotides for encoding.
    :param segment_length: Segments should be at most this length.
    :return: Encoded sequence as integers.
    """
    encoded_seqs = []
    for sequence in sequences:
        encoded_seqs.append(encode_sequence(sequence, nuc_order))

    #return np.array(pad_sequences(encoded_seqs, maxlen=segment_length, padding="post"))
    return np.array(encoded_seqs)

def load_virus_model(model_path):
    model = load_model(model_path)
    return model


def proc_fastq(infile, target_length):
    """
    Parse a FASTQ file and extract reads. Rejects reads that:
        1) Are an unexpected length
        2) Contain > 1 N
        3) Contain a character other than ACGTN
    :param infile: path to a FASTQ file.
    :param target_length: expected length of reads. Reads of other lengths are rejected
    :return: encoded_c, Accepted reads prepared for model use as follows:
        1) N's replacted with ranom nucleotide
        2) Padded (randomly) or truncated to 76bp
        Represent with integers 0-3
    :return: raw_seqs, Accepted reads without modification, for assembly
    """
    #Determine padding to acieve 76bp (SEGMENT_LENGTH) length
    length_difference = SEGMENT_LENGTH - target_length
    end_ns = "N" * length_difference
    #Validity REGEXs
    #note: triple braces: outer two are converted to single leteral for regex, inner evaluates variable
    seq_remove = re.compile(f"(?:[AN]{{{SEARCHSUBLEN}}})|(?:[CN]{{{SEARCHSUBLEN}}})|(?:[GN]{{{SEARCHSUBLEN}}})|(?:[TN]{{{SEARCHSUBLEN}}})")
    seq_valid = re.compile(f"[ACTGN]{{{target_length}}}")
    
    #Parse FASTQ
    f = open(infile)
    l = f.readlines()
    f.close()
    all_seqs = [l[i].strip() for i in range(1, len(l), 4)]
    #print("Total reads read:", len(all_seqs))
    seqs = list(np.unique(all_seqs))
    #print("Unique reads read:", len(seqs))
    raw_seqs = []
    processed_seqs = []
    #random.seed(0)
    for seq in seqs:
        #THIS CODE IS FOR RANDOMLY ADDING N's IN TESTING
        #if random.random() < 0.05:
        #    idx = random.randrange(target_length)
        #    seq = seq[:idx] + "N" + seq[idx+1:]
        #END TEST CODE    
        if re.fullmatch(seq_valid, seq) and seq.count("N") < 2 and re.search(seq_remove, seq) is None:
            raw_seqs.append(seq)
            if length_difference > 0:
                seq += end_ns
            processed_seqs.append(handle_N(seq[:SEGMENT_LENGTH]))
                      
    #Now done earler, on raw sequences - a bit inefficient, and counts Ns as letters
    #seqs = list(np.unique([handle_non_ATGC(i) for i in seqs]))
    #Replaced by a single fixed size required of all reads, else they are excluded
    #medsize = np.median([len(i) for i in seqs])
    #seqs = [i[:int(medsize)] for i in seqs]
    
    encoded_c = encode_sequences(processed_seqs)

    return encoded_c, raw_seqs

def make_clist(lst):
    return (c_char_p * len(lst))(*[x.encode() for x in lst])

def assemble_read_call_c(bact_scores, viral_scores, reads, seed_indices, filen):
    """
    ToDo: Daniel, please document better
    :param bact_scores: all reads assigned bact scores
    :param viral_scores: all reads assigned viral scores
    :param reads: all reads
    :param seed_indices: the bacterial and viral seeds with score>0.7 for either
    :param filen: file name
    :return:
    """
    librd = CDLL("/wistar/auslander/Daniel/scripts/src/assemble_read_c.so")
    librd.connect()

    # preparing these variables for c
    arr_f1 = (c_float * len(bact_scores))(*list(bact_scores))
    arr_f2 = (c_float * len(viral_scores))(*list(viral_scores))
    arr_s = (c_int * len(seed_indices))(*list(seed_indices))
    filen = c_char_p(filen.encode())
    arr_ch = make_clist(reads)
    m = c_int(len(reads))
    s = c_int(len(seed_indices))
    n = c_int(len(reads[0])) #Maybe should be a constant?
    #librd.assemble_read.restype = c_char_p
    #print("c called")
    #result = librd.assemble_read_loop(arr_f1, arr_f2, arr_ch, arr_s, n, m, s, filen)
    #print(datetime.now().time())
    #yep.start("output.prof")

    #Passed to C: bact & viral scores, all reads, seed indices, len(read), # of reads, # of seeds, output file
    print('calling C', flush=True)
    librd.assemble_read_loop(arr_f1, arr_f2, arr_ch, arr_s, n, m, s, filen)
    #yep.stop()
    #print(datetime.now().time())
    #return result



def extract_contigs(invars):
    '''

    :param invars: input file, output file, read length, model path
    :return:
    '''

    inpath = invars[0]
    outpath = invars[1]
    target_length = invars[2]
    model_path = invars[3]

    # output file name
    fn = outpath + inpath.split('/')[-1].replace('_unmapped', '_contigs').replace('fastq', 'txt')
    if '_contigs' not in fn:
        fn = fn.replace('.txt', '_contigs.txt')

    file_exists = os.path.exists(fn)
    if file_exists:
        return 0
    else:

        ##Thats the trained neural network model
        model = load_virus_model(model_path)
        #print("loaded model")
        #Encode unmapped RNAseq reads for model
        encoded_seqs, raw_seqs = proc_fastq(inpath, target_length)
        #print("processed FASTQ. found", len(raw_seqs), "seqs")
        # Predict the viral unmapped RNAseq reads using the model
        if len(raw_seqs) > 0:
            scc = model.predict(encoded_seqs)
            #print("ran model")
            # Select seeds for assembling contigs - reads scored more than 0.7 by model
            qq = 0.46
            #Just a naive implementation of conditions
            seed_bact_score_indices = [(scc[i,2],i) for i in range(scc.shape[0]) if scc[i,2] > max(qq, scc[i,1], scc[i,0])]
            seed_virus_score_indices = [(scc[i,0],i) for i in range(scc.shape[0]) if scc[i,0] >  max(qq, scc[i,1], scc[i,2])]
            seed_bact_score_indices.sort(key=lambda x : x[0], reverse=True)
            seed_virus_score_indices.sort(key=lambda x : x[0], reverse=True)
            
            #bact_sorted_indices = np.argsort(scc[:,2])[::-1] ##sorted indices from greater to smaller by bacteria score
            #conf_bact_indices = bact_sorted_indices[:np.count_nonzero(scc[:,2] > qq)] #This is indices of bac>0.7 sorted by bacterial score
            #virus_sorted_indices = np.argsort(scc[:,0])[::-1] ##sorted indices from greater to smaller by virus score
            #conf_virus_indices = virus_sorted_indices[:np.count_nonzero(scc[:,0] > qq)]#This is indices of vir>0.7 sorted by bacterial score

            #print("ci", len(conf_bact_indices), len(conf_virus_indices))
            #print("calling c")
            #print(list(bact_sorted_indices) + list(virus_sorted_indices))
            assemble_read_call_c(list(scc[:,2]), list(scc[:,0]), raw_seqs, [x[1] for x in seed_bact_score_indices] + [x[1] for x in seed_virus_score_indices], fn)
        else:
            print("Found no valid sequences in", inpath)
        return 0
        

def run_virna_pred(inpath,outpath,multi_proc,model_path,num_threads, target_length=76):


    ##get input fastq files
    infastq = list(set([i for i in glob.glob(inpath + '*.fastq') if '.fastq' in i]))

    ##To make sure the outputs were not yet generated, get the input fastq that do not have output contigs yet
    outs = glob.glob(outpath + '/' + '*.txt')
    nmi = [i.split('/')[-1].replace('.fastq', '').replace('_unmapped', '') for i in infastq]
    nmo = [i.split('/')[-1].replace('.txt', '').replace('_contigs', '') for i in outs]
    infiles = [infastq[i] for i in range(len(nmi)) if nmi[i] not in nmo]
    print('found', len(infiles), 'input files')
    print('starting_prediction...', flush=True)
    if multi_proc:
        freeze_support()
        pool = Pool(processes=num_threads)
        pool.map(extract_contigs, [[f, outpath, target_length, model_path] for f in infiles])

    else:
        ins = [[f, outpath, target_length, model_path] for f in infiles]
        for i in range(len(ins)):
            extract_contigs(ins[i])

    #print("Done processing")



def main():
    #print(datetime.now().time())
    PWD = os.getcwd()
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--input", type=str, help="input directory", required=True
    )
    parser.add_argument(
        "--output", type=str, help="output directory", required=True
    )
    parser.add_argument(
        "--read_length", type=int, help="Expected length of RNAseq reads", required=False, default=76
    )
    parser.add_argument("--num_threads", type=int, help="number of threads to run with pool multi processing", required=False, default=1
    )
    parser.add_argument(
        "--model_path", type=str, help="path to Tensorflow model to predict whether reads of a fixed length come from viruses or not", required=True
    )

    args = parser.parse_args()
    print("a", flush=True)
    inpath = args.input
    outpath = args.output

    if not os.path.isdir(inpath):
        print('input directory %s not found',inpath)
        sys.exit(1)
    #print("b", flush=True)
    if not os.path.isdir(outpath):
        print('output directory %s not found',outpath)
        sys.exit(1)
    #print("c", flush=True)
    if args.num_threads == 1:
        multi_proc = False
        num_threads = 1
    elif args.num_threads > 1:
        num_threads = args.num_threads
        multi_proc  = True
    #print("d", flush=True)
    if not os.path.isdir(args.model_path):
        print(f'model {args.model_path} not found')
        sys.exit(1)
    else:
        model_path = args.model_path
    #print("e", flush=True)

    print("Reading fastq at {}...".format(inpath), flush=True)
    run_virna_pred(inpath, outpath, multi_proc, model_path,num_threads, args.read_length)
    #print("done")
    #print(datetime.now().time())

if __name__ == '__main__':
    #print("z", flush=True)
    main()
