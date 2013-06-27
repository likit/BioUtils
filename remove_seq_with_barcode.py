'''Remove sequences containing a barcode in input_barcode

input_barcode can be a barcode sequence or a file containing a list of
barcodes.

The script assume input sequences in FASTA format with one line per sequence.

'''


import sys
import os


def read_barcode(input_barcode, is_file=False):
    barcodes = set()
    if is_file:
        for line in open(input_barcode):
            barcodes.add(line.strip())
    else:
        if input_barcode:
            barcodes.add(input_barcode)
            print >> sys.stderr, 'barcodes = %s' % input_barcode

    return barcodes


def remove(filename, barcodes):
    fp = open(filename)
    while (1):
        try:
            line = fp.next()
        except StopIteration:
            break;
        else:
            if line.startswith('>'):
                newbc = line.split()[3].split('=')[1]
                if newbc in barcodes:
                    line = fp.next()
                    continue
                else:
                    print line.strip()  # seqeunce ID
                    line = fp.next()
                    print line.strip()  # sequence


def main():
    input_seq = sys.argv[1]
    input_barcode = sys.argv[2]
    barcodes = read_barcode(input_barcode,
            is_file=os.path.isfile(input_barcode))
    remove(input_seq, barcodes)


if __name__=='__main__':
    main()
