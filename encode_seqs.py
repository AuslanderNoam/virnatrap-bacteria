
import sys
import numpy as np
from tqdm import tqdm
from Bio import SeqIO


#Outfile, then one or more FASTA



def onehot_seq(seq):
    letter_to_index =  {'A':0, 'a':0,
                        'C':1, 'c':1,
                        'G':2, 'g':2,
                        'T':3, 't':3}
    to_return = np.zeros((len(seq),4), dtype='int8')
    for idx,letter in enumerate(seq):
      if letter in letter_to_index:
        to_return[idx,letter_to_index[letter]] = 1
    return to_return

def encode_sequence(fasta, shuffle = False):
    x = np.array([onehot_seq(seq) for seq in tqdm(SeqIO.parse(fasta, "fasta"))])
    #No reverse complement
    print(f'There are {x.shape[0]} examples in this class')
    #xcon = np.concatenate((x_pos, x_neg))
    #This adds a fourth dimension so each elem is in its own array
    #It's needed if will be using a 2D conv. 
    #x = np.expand_dims(xcon, axis=3)
    if shuffle:
        indices = np.arange(x.shape[0])
        np.random.shuffle(indices)
        x = x[indices,:]
    return x

i = 0
xs = []
count = 0
precounts = []
for name in sys.argv[2:]:
    print("loading file", name, "with label", i)
    x = encode_sequence(name, shuffle = False)
    xs.append(x)
    count += x.shape[0]
    precounts.append(count)
    i += 1

x_full = np.concatenate(xs)
y_full = np.zeros(count)
#Assign consecutive labels
for pc in precounts:
    y_full[pc:] += 1

np.savez_compressed(sys.argv[1], x=x_full, y=y_full, allow_pickle=True)

    
