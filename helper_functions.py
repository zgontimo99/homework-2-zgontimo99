def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("the brown cat", "these brownies", lambda x, y: [-1, 1][x == y])
    ('the-- brown cat', 'these brownies-', 3.0)

    Other alignments are also possible.

    """
    n = len(seq1)
    m = len(seq2)
    dp = [[0] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        dp[i][0] = dp[i - 1][0] + scoring_function(seq1[i - 1],'*')
    for j in range(1, m + 1):
        dp[0][j] = dp[0][j - 1] + scoring_function('*',seq2[j - 1])

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            if seq1[i - 1] == seq2[j - 1]:
                match_score = dp[i - 1][j - 1] + scoring_function(seq1[i - 1],seq2[j - 1])
            else:
                match_score = dp[i - 1][j - 1] + scoring_function(seq1[i - 1],seq2[j - 1])
            gap1_score = dp[i - 1][j] + scoring_function(seq1[i - 1],'*')
            gap2_score = dp[i][j - 1] + scoring_function('*',seq2[j - 1])
            dp[i][j] = max(match_score, gap1_score, gap2_score)

    
    align1, align2 = '', ''
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and dp[i][j] == dp[i - 1][j] + scoring_function(seq1[i - 1],'*'):
            align1 = seq1[i - 1] + align1
            align2 = '-' + align2
            i -= 1
        elif j > 0 and dp[i][j] == dp[i][j - 1] + scoring_function('*',seq2[j - 1]):
            align1 = '-' + align1
            align2 = seq2[j - 1] + align2
            j -= 1
        else:
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1

    return align1, align2, dp[n][m]


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("the brown cat", "these brownies", lambda x, y: [-1, 1][x == y])
    ('the-- brown', 'these brown', 7.0)

    Other alignments are also possible.

    """
    n = len(seq1)
    m = len(seq2)
    score_matrix = [[0] * (m + 1) for _ in range(n + 1)]
    traceback_matrix = [[0] * (m + 1) for _ in range(n + 1)]

    max_score = 0
    max_i, max_j = 0, 0

    for i in range(1, n+1):
        for j in range(1, m+1):
            match = score_matrix[i - 1][j - 1] + scoring_function(seq1[i - 1],seq2[j - 1])
            delete = score_matrix[i - 1][j] + scoring_function(seq1[i - 1],'*')
            insert = score_matrix[i][j - 1] + scoring_function('*',seq2[j - 1])
            score_matrix[i][j] = max(0, match, delete, insert)

            if score_matrix[i][j] == 0:
                traceback_matrix[i][j] = 0
            elif score_matrix[i][j] == match:
                traceback_matrix[i][j] = 1  # Diagonal
            elif score_matrix[i][j] == delete:
                traceback_matrix[i][j] = 2  # Up
            elif score_matrix[i][j] == insert:
                traceback_matrix[i][j] = 3  # Left

            if score_matrix[i][j] > max_score:
                max_score = score_matrix[i][j]
                max_i, max_j = i, j

    
    aligned_sequence1 = []
    aligned_sequence2 = []
    i, j = max_i, max_j
    while traceback_matrix[i][j] != 0:
        if traceback_matrix[i][j] == 1:  # Diagonal
            aligned_sequence1.append(seq1[i - 1])
            aligned_sequence2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == 2:  # Up
            aligned_sequence1.append(seq1[i - 1])
            aligned_sequence2.append('-')
            i -= 1
        elif traceback_matrix[i][j] == 3:  # Left
            aligned_sequence1.append('-')
            aligned_sequence2.append(seq2[j - 1])
            j -= 1

    aligned_sequence1.reverse()
    aligned_sequence2.reverse()

    return ''.join(aligned_sequence1), ''.join(aligned_sequence2), max_score
