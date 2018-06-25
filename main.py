import sys
import alignement
import numpy
import collections

seq1 = 'MGEIGFTEKQEALVKESWEILKQDIPKYSLHFFSQILEIAPAAKGLFSFLRDSDEVPHNNPKLKAHAVKVFKMTCETAIQLREEGKVVVADTTLQYLGSIHLKSGVIDPHFEVVKEALLRTLKEGGEKYNEEVEGAWSQAYDHLALAIKTEMKQEES'

seq2 = 'MEKVPGEMEIERRERSEELSEAERKAVQATWARLYANCEDVGVAILVRFFVNFPSAKQYFSQFKHMEEPLEMERSPQLRKHACRVMGALNTVVENLHDPEKVSSVLSLVGKAHALKHKVEPVYFKLSGVILEVIAEEFANDFPPETQRAWAKLRGLIYSHVTAAYKEVGWVQQVPNATTPPATLPSSGP'

seq3 = 'MVLSAADKNNVKGIFTKIAGHAEEYGAETLERMFTTYPPTKTYFPHFDLSHGSAQIKGHGKKVVAALIEAANHIDDIAGTLSKLSDLHAHKLRVDPVNFKLLGQCFLVVVAIHHPAALTPEVHASLDKFLCAVGTVLTAKYR'

seq4 = 'MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKHPGDFGADAQGAMNKALELFRKDMASNYKELGFQG'

seq5 = 'MERLESELIRQSWRAVSRSPLEHGTVLFSRLFALEPSLLPLFQYNGRQFSSPEDCLSSPEFLDHIRKVMLVIDAAVTNVEDLSSLEEYLATLGRKHRAVGVRLSSFSTVGESLLYMLEKCLGPDFTPATRTAWSQLYGAVVQAMSRGWDGE'

def calcul_sp():
    sp = 0
    with open('alignement.txt') as f:
        sequenceFileContent = f.readlines()
    print sequenceFileContent
    for i in range (len(sequenceFileContent)):
        for j in range (len(sequenceFileContent)):
            for x in range (len(sequenceFileContent[i])):
                if (sequenceFileContent[i][x] != sequenceFileContent[j][x]):
                    sp = sp + 1
    return sp

def sequence_consensus():
    with open('alignement.txt') as f:
        sequenceFileContent = f.readlines()
    totalstr = ''
    for i in range (len(sequenceFileContent[0])):
        str = sequenceFileContent[0][i] + sequenceFileContent[1][i] + sequenceFileContent[2][i] + sequenceFileContent[3][i] + sequenceFileContent[4][i]
        totalstr = totalstr + collections.Counter(str).most_common(1)[0][0]
        str = ''
    print totalstr

sequence_consensus()

alignement.read_blosum()

if (len(sys.argv) == 3):
    print("Sequence 1 : " + sys.argv[1])
    print("Sequence 2 : " + sys.argv[2])
    [E, F, G] = alignement.matrice_distance(sys.argv[1], sys.argv[2])
    [str1, str2, score] = alignement.backtracking(sys.argv[1], sys.argv[2], E, F, G)
    print("Alignement obtenu :")
    print(str1)
    print(str2)
    print("\n")
else:
    try:
        with open('sequences.fasta') as f:
            sequenceFileContent = f.readlines()
    except:
        print 'Problem with file. Is the file sequences.fasta in the same folder as this python file ?'
    sequences = []
    sequenceNumber = 0
    seq = ''
    for i in range(len(sequenceFileContent)):
        if (sequenceFileContent[i][0] == '>'):
            if (sequenceNumber != 0):
                sequences.append(seq)
            sequenceNumber = sequenceNumber + 1
            seq = ''
        else:
            seq = seq + sequenceFileContent[i].replace('\n', '')
    if (seq != ''):
        sequences.append(seq)
    scoreEtoile = numpy.zeros((len(sequences), len(sequences))).astype(int)
    for i in range(len(sequences)):
        for j in range(len(sequences)):
            print("Sequence 1 : " + sequences[i])
            print("Sequence 2 : " + sequences[j])
            [E, F, G] = alignement.matrice_distance(sequences[i], sequences[j])
            [str1, str2, score] = alignement.backtracking(sequences[i], sequences[j], E, F, G)
            print("Alignement obtenu :")
            print(str1)
            print(str2)
            print("\n")
            scoreEtoile[i][j] = score
    print scoreEtoile
