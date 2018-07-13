import sys
import random


def int_to_base(number, fA, fC, fG, fT):
    if number < fA:
        return 'A'
    elif number < fA + fC:
        return 'C'
    elif number < fA + fC + fG:
        return 'G'
    elif number < fA + fC + fG + fT:
        return 'T'
    else:
        return 'N'


def readGenerator(fA, fC, fG, fT, n_reads, read_length, output_filename, seq_id=""):
    if seq_id is "":
        read_id = ">RANDOMREAD."

    r_len = int(read_length)
    percentage = int(n_reads)/100
    read_seq = ""

    with open(output_filename, 'w') as f:
        for i in range(1, int(n_reads)+1, 1):
            f.write(read_id + str(i) + "\n")
            for j in range(1,r_len+1,1):
                next_base = int_to_base(random.randint(0, 99), int(fA), int(fC), int(fG), int(fT))
                read_seq += next_base
                if j % 70 == 0:
                    f.write(read_seq + '\n')
                    read_seq = ""
            f.write(read_seq + '\n')
            read_seq = ""
            if i % percentage == 0:
                print("Iter: " + str(i*100/int(n_reads)))

# ------------------------

def main():
    if len(sys.argv) < 8 or len(sys.argv) > 9 :
        #  python randomReadsGenerator.py 25 25 25 25 100 130 random_test.fa
        #                                          1    2    3    4         5                   6                  7           8
        print("USAGE: python randomReadsGenerator.py <fA> <fC> <fG> <fT> <Number of Reads> <Read_Average_Length> <output> <OPTIONAL-id>")
        exit(-1)

    if len(sys.argv) == 8:
        readGenerator(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7])
    elif len(sys.argv) == 9:
        readGenerator(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], sys.argv[7], sys.argv[8])

# ----------------------


if __name__ == "__main__":
    main()
