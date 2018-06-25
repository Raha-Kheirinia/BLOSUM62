global match
global mismatch
global MIN
global A
global B

score = 0
A   = -10.
B   = -1.
match = 1.
mismatch = -0.
MIN = -float("inf")

def read_blosum():
    with open('BLOSUM62.txt') as f:
        fileContent = f.readlines()
    global alignment_score
    alignment_score = []
    for i in range (1, len(fileContent)):
        alignment_score = alignment_score + [ fileContent[i].split() ]
    alignment_score = [line[1:] for line in alignment_score]
    global alignmentIndices
    alignmentIndices = fileContent[0].split()

def get_score(a, b):
    indice_a = -1
    indice_b = -1
    for i in range(0, len(alignmentIndices)):
        if (a == alignmentIndices[i]):
            indice_a = i
        if (b == alignmentIndices[i]):
            indice_b = i
    return int(alignment_score[indice_a][indice_b])

def backtracking(s, t, E, F, G):
    seq1 = ''
    seq2 = ''
    i = len(t)
    j = len(s)
    score = G[i][j]
    print score
    while (i>0 or j>0):
        if (i>0 and j>0 and G[i][j] == G[i-1][j-1] + get_score(t[i - 1], s[j - 1])):
            seq1 += s[j-1]
            seq2 += t[i-1]
            i -= 1; j -= 1
        elif (i>0 and G[i][j] == F[i][j]):
            seq1 += '_'
            seq2 += t[i-1]
            i -= 1
        elif (j>0 and G[i][j] == E[i][j]):
            seq1 += s[j-1]
            seq2 += '_'
            j -= 1

    seq1r = ' '.join([seq1[j] for j in range(-1, -(len(seq1)+1), -1)])
    seq2r = ' '.join([seq2[j] for j in range(-1, -(len(seq2)+1), -1)])

    return [seq1r, seq2r, score]

def matrice_distance(s, t):
    E = [[initialisation_x(i, j) for j in range(0, len(s) + 1)] for i in range(0, len(t) + 1)]
    F = [[initialisation_y(i, j) for j in range(0, len(s) + 1)] for i in range(0, len(t) + 1)]
    G = [[initialisation_m(i, j) for j in range(0, len(s) + 1)] for i in range(0, len(t) + 1)]
    for j in range(1, len(s) + 1):
        for i in range(1, len(t) + 1):
            E[i][j] = max((A + B + G[i][j-1]), (B + E[i][j-1]), (A + B + F[i][j-1]))
            F[i][j] = max((A + B + G[i-1][j]), (A + B + E[i-1][j]), (B + F[i-1][j]))
            G[i][j] = max(get_score(t[i - 1], s[j - 1]) + G[i-1][j-1], E[i][j], F[i][j])
    return [E, F, G]

def initialisation_x(i, j):
    if i > 0 and j == 0:
        return MIN
    else:
        if j > 0:
            return -10 + (-0.5 * j)
        else:
            return 0

def initialisation_y(i, j):
    if j > 0 and i == 0:
        return MIN
    else:
        if i > 0:
            return -10 + (-0.5 * i)
        else:
            return 0

def initialisation_m(i, j):
    if j == 0 and i == 0:
        return 0
    else:
        if j == 0 or i == 0:
            return MIN
        else:
            return 0
