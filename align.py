class GlobalSeqAlignment:
    def __init__(self, x: str, y: str, M=4, m=-2, g=-2):
        self.x, self.y = x, y
        self.M, self.m, self.g = M, m, g
        self.score = []
        for i in range(len(x) + 1):
            row = [0] * (len(y) + 1)
            self.score.append(row)
        self.traceback = []
        for i in range(len(x) + 1):
            row = [''] * (len(y) + 1)
            self.traceback.append(row)
        self.initialize_matrices()
        self.calculate_matrices()
        self.print_score_matrix()
        self.print_traceback_matrix()

    def initialize_matrices(self):
        for j in range(1, len(self.y) + 1):
            self.score[0][j] = self.g * j
        for i in range(1, len(self.x) + 1):
            self.score[i][0] = self.g * i

    def calculate_matrices(self):
        for i in range(1, len(self.x) + 1):
            for j in range(1, len(self.y) + 1):
                left = self.score[i][j - 1] + self.g
                up = self.score[i - 1][j] + self.g
                if self.x[i - 1] == self.y[j - 1]:
                    match_score = self.M
                else:
                    match_score = self.m
                diag = self.score[i - 1][j - 1] + match_score

                max_score = max(left, up, diag)
                self.score[i][j] = max_score

                if max_score == diag:
                    self.traceback[i][j] = 'D'
                elif max_score == left:
                    self.traceback[i][j] = '<'
                else:
                    self.traceback[i][j] = '^'

                if max_score == diag:
                    if left == diag or up == diag:
                        self.traceback[i][j] = 'D'
                    elif left > diag or up > diag:
                        self.traceback[i][j] = '<' if left > up else '^'
                elif max_score == left:
                    if up == left:
                        self.traceback[i][j] = '^'

    def print_score_matrix(self):
        print("Score Matrix:")
        max_width = max(len(str(max(row)) + " " + str(min(row))) for row in self.score)
        for row in self.score:
            formatted_row = [f"{num:>{max_width}}" for num in row]
            print(" ".join(formatted_row))

    def print_traceback_matrix(self):
        print("Traceback Matrix:")
        for row in self.traceback:
            print(" ".join(row))

    def get_optimal_alignment(self):
        i, j = len(self.x), len(self.y)
        x_aligned, y_aligned = '', ''

        while i > 0 or j > 0:
            if self.traceback[i][j] == 'D':
                x_aligned = self.x[i - 1] + x_aligned
                y_aligned = self.y[j - 1] + y_aligned
                i -= 1
                j -= 1
            elif self.traceback[i][j] == '<':
                x_aligned = '-' + x_aligned
                y_aligned = self.y[j - 1] + y_aligned
                j -= 1
            else:
                x_aligned = self.x[i - 1] + x_aligned
                y_aligned = '-' + y_aligned
                i -= 1

        return x_aligned, y_aligned


def parse(filename):
    sequences = []
    sequence = ''
    with open(filename, 'r') as fasta_file:
        for line in fasta_file:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    sequences.append(sequence)
                    sequence = ''
                if len(sequences) == 2:
                    break
            else:
                sequence += line
        if sequence and len(sequences) < 2:
            sequences.append(sequence)
    while len(sequences) < 2:
        sequences.append('')
    return sequences[0], sequences[1]


def parse_FASTA_file(fasta_stream):
    sequences = []
    sequence = ''
    for line in fasta_stream:
        line = line.strip()
        if line.startswith(">"):
            if sequence:
                sequences.append(sequence)
                sequence = ''
            if len(sequences) == 2:
                break
        else:
            sequence += line
    if sequence and len(sequences) < 2:
        sequences.append(sequence)
    while len(sequences) < 2:
        sequences.append('')
    return sequences[0], sequences[1]


# def main():
#     x, y = parse_FASTA_file(sys.stdin)
#     align = GlobalSeqAlignment(x, y)
#     x_prime, y_prime = align.get_optimal_alignment()
#     print(x_prime)
#     print(y_prime)
def write_alignment_to_file(x_prime, y_prime, filename="my.output1"):
    with open(filename, 'w') as file:
        file.write(x_prime + '\n')
        file.write(y_prime + '\n')


def main():
    x, y = parse("aligntest.input2.txt")  # Assuming parse is a function that reads and returns the sequences
    align = GlobalSeqAlignment(x, y)  # Assuming GlobalSeqAlignment is initialized properly
    x_prime, y_prime = align.get_optimal_alignment()  # Assuming this method returns the aligned sequences
    print(x_prime)
    print(y_prime)
    write_alignment_to_file(x_prime, y_prime)


if __name__ == "__main__":
    main()
