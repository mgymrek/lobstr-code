import sys

match=1.
mismatch=-1
gap=-2
s={
'AA':     match,'AG':mismatch,'AC':mismatch,'AT':mismatch,\
'GA':mismatch,'GG':     match,'GC':mismatch,'GT':mismatch,\
'CA':mismatch,'CG':mismatch,'CC':     match,'CT':mismatch,\
'TA':mismatch,'TG':mismatch,'TC':mismatch,'TT':     match,\
}

def nw_noeg(seq1, seq2, aln=False):
    """
    Global alignment of two sequences with
    no end-gap penalty.
    Normalize by length of overlap
    Returns alignment score
    """
    rows=len(seq1)+1
    cols=len(seq2)+1
    # init matrix
    a=[]
    for i in range(rows):
        a+=[[0.]*cols]

    for i in range(rows):
        a[i][0] = 0
    for j in range(cols):
        a[0][j] = 0

    # Run DP
    for i in range(1,rows):
        for j in range(1,cols):
            choice1 = a[i-1][j-1] + s.get((seq1[i-1] + seq2[j-1]),mismatch)
            choice2 = a[i-1][j] + gap
            choice3 = a[i][j-1] + gap
            a[i][j] = max(choice1, choice2, choice3)

    # calculate nucs overlapping
    maxscore = 0
    maxoverlap = 0
    maxbottomrow = max(a[i])
    maxrightcol = max([a[k][j] for k in range(rows)])
    if maxbottomrow > maxrightcol:
        maxscore = maxbottomrow
        bottom_index = a[i].index(maxscore)
        # ends at last position of seq 1 and position bottom_index of seq2
        maxoverlap = min(i, bottom_index)
    else:
        maxscore = maxrightcol
        right_index = [a[k][j] for k in range(rows)].index(maxscore)
        maxoverlap = min(j, right_index)
    try:
        score = maxscore*1.0/maxoverlap
    except: score = 0
    if not aln: return score

    # get aligned sequences				
    aseq1 = ''
    aseq2 = ''
    i = len(seq1)
    j = len(seq2)
    while i>0 and j>0:
        score = a[i][j]
        score_diag = a[i-1][j-1]
        score_up = a[i][j-1]
        score_left = a[i-1][j]
        if score == score_diag + s[seq1[i-1] + seq2[j-1]]:
            aseq1 = seq1[i-1] + aseq1
            aseq2 = seq2[j-1] + aseq2
            i -= 1
            j -= 1
        elif score == score_left + gap:
            aseq1 = seq1[i-1] + aseq1
            aseq2 = '-' + aseq2
            i -= 1
        elif score == score_up + gap:
            aseq1 = '-' + aseq1
            aseq2 = seq2[j-1] + aseq2
            j -= 1
        else:
            #should never get here..
            i=0
            j=0
            aseq1='ERROR';aseq2='ERROR';seq1='ERROR';seq2='ERROR'
    while i>0:
    #If we hit j==0 before i==0 we keep going in i.
        aseq1 = seq1[i-1] + aseq1
        aseq2 = '-' + aseq2
        i -= 1		

    while j>0:
    #If we hit i==0 before i==0 we keep going in j. 
        aseq1 = '-' + aseq1
        aseq2 = seq2[j-1] + aseq2
        j -= 1
    return score, aseq1, aseq2
