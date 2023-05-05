
import sys, random
from Bio import SeqIO


def seq2reads(item, size, step):
    reads = []
    seq = item.seq
    l = len(seq)
    for i in range(0, l, step):
        j = i + size
        if j <= l:
            text = str(str(seq[i:j]))
            if "N" in text:
                print("excluded read with N", file=sys.stderr)
            else:
                reads.append(">" + ";".join((item.id, str(i+1), str(j))) + "\n" + text + "\n")
    return reads

size = int(sys.argv[2])
step = int(sys.argv[3])
target_train = int(sys.argv[4])
target_val = int(sys.argv[5])

seqs = SeqIO.parse(sys.argv[1], "fasta")
items = []
for item in seqs:
    items.append(item)

#seed = random.randrange(sys.maxsize)
#rng = random.Random(seed) #NOTE: before, claimed seed was not actually used.
#random.seed(seed)
#print("Seed was:", seed)
random.shuffle(items)
print(len(items))
train_count = 0
val_count = 0
test_count = 0
state = 0

trainFile = open(sys.argv[6] + "_train.txt", "w")
valFile = open(sys.argv[6] + "_val.txt", "w")
testFile = open(sys.argv[6] + "_test.txt", "w")

for item in items:
    reads = seq2reads(item, size, step)
    if state == 0:
        trainFile.write("".join(reads))
        train_count += len(reads)
        if train_count > target_train:
            state = 1
    elif state == 1:
        valFile.write("".join(reads))
        val_count += len(reads)
        if val_count > target_val:
            state = 2
    elif state == 2:
        testFile.write("".join(reads))
        test_count += len(reads)
trainFile.close()
valFile.close()
testFile.close()
print(train_count)
print(val_count)
print(test_count)    
    
    
    
    
    
    
    
    
    
    
    


    
