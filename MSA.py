import argparse
import re
import numpy as np
from Bio import SeqIO
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib import colors
import tkinter as tk
from tkinter import filedialog
from tkinter import messagebox


def is_valid_dna(seq):
    """
    Validates that the sequence contains only ACGT characters.

    Returns:
        bool: If the sequence contains only ACGT characters.
    """
    return bool(re.fullmatch(r"[ACGTacgt]+", seq))


def parse_sequences(file=None, seqs=None):
    """
    Parses sequences from either file or direct input.

    Returns:
        sequences (list): A list of sequences.
    """
    sequences = []
    if file:
        try:
            for record in SeqIO.parse(file, "fasta"):
                seq = str(record.seq).upper()
                if is_valid_dna(seq):
                    sequences.append(seq)
            if len(sequences) < 2:
                raise ValueError("FASTA file must contain at least two valid DNA sequences.")
        except Exception as e:
            raise ValueError(f"Error reading FASTA file: {e}")

    elif seqs:
        for seq in seqs:
            if is_valid_dna(seq):
                sequences.append(seq.upper())
            else:
                raise ValueError(f"Invalid DNA sequence provided: {seq}")
        if len(sequences) < 2:
            raise ValueError("At least two valid sequences must be provided via --seqs.")

    else:
        raise ValueError("You must provide either --file or --seqs input.")

    return sequences


def initialize_matrix(len_a, len_b, gap_penalty):
    """
    Initializes the matrix of size (len_a, len_b).

    Arguments:
        len_a (int): The length of the first sequence.
        len_b (int): The length of the second sequence.
        gap_penalty (int): The gap penalty.

    Returns:
        matrix (array): The matrix of size (len_a, len_b).
    """
    matrix = np.zeros((len_a + 1, len_b + 1))
    for i in range(1, len_a + 1):
        matrix[i][0] = matrix[i - 1][0] + gap_penalty
    for j in range(1, len_b + 1):
        matrix[0][j] = matrix[0][j - 1] + gap_penalty
    return matrix


def compute_scores(a, b, matrix, match_score, mismatch_score, gap_penalty):
    """
    Computes the scores between two sequences.

    Arguments:
        a (int): The first sequence.
        b (int): The second sequence.
        matrix (array): The matrix of size (len_a, len_b).
        match_score (int): The match score.
        mismatch_score (int): The mismatch score.
        gap_penalty (int): The gap penalty.
    """
    for i in range(1, len(a) + 1):
        for j in range(1, len(b) + 1):
            score = match_score if a[i - 1] == b[j - 1] else mismatch_score
            matrix[i][j] = max(
                matrix[i - 1][j] + gap_penalty,
                matrix[i][j - 1] + gap_penalty,
                matrix[i - 1][j - 1] + score
            )


def traceback(a, b, matrix, match_score, mismatch_score, gap_penalty):
    """
    Computes the traceback between two sequences.

    Arguments:
        a (int): The first sequence.
        b (int): The second sequence.
        matrix (array): The matrix of size (len_a, len_b).
        match_score (int): The match score.
        mismatch_score (int): The mismatch score.
        gap_penalty (int): The gap penalty.

    Returns:
        aligned_a,aligned_b (int): The aligned sequences.
    """
    aligned_a, aligned_b = "", ""
    i, j = len(a), len(b)
    while i > 0 or j > 0:
        if i > 0 and j > 0 and matrix[i][j] == matrix[i - 1][j - 1] + (
            match_score if a[i - 1] == b[j - 1] else mismatch_score):
            aligned_a = a[i - 1] + aligned_a
            aligned_b = b[j - 1] + aligned_b
            i -= 1
            j -= 1
        elif i > 0 and matrix[i][j] == matrix[i - 1][j] + gap_penalty:
            aligned_a = a[i - 1] + aligned_a
            aligned_b = "-" + aligned_b
            i -= 1
        else:
            aligned_a = "-" + aligned_a
            aligned_b = b[j - 1] + aligned_b
            j -= 1
    return aligned_a, aligned_b


def needleman_wunsch(a, b, match_score=1, mismatch_score=0, gap_penalty=-1):
    """
    Performs the Needleman-Wunsch between two sequences.

    Arguments:
        a (int): The first sequence.
        b (int): The second sequence.
        match_score (int): The match score.
        mismatch_score (int): The mismatch score.
        gap_penalty (int): The gap penalty.

    Returns:
        aligned_a,aligned_b (int): The aligned sequences.
    """
    matrix = initialize_matrix(len(a), len(b), gap_penalty)
    compute_scores(a, b, matrix, match_score, mismatch_score, gap_penalty)
    aligned_a, aligned_b = traceback(a, b, matrix, match_score, mismatch_score, gap_penalty)
    return aligned_a, aligned_b


def center_star_msa(sequences, match_score=1, mismatch_score=0, gap_penalty=-1):
    """
    Performs Center-Star multiple sequence alignment on a list of DNA sequences.

    Args:
        sequences (list): A list of DNA sequences to align.
        match_score (int): The score for a match between two bases. Default is 1.
        mismatch_score (int): The penalty for a mismatch between two bases. Default is -1.
        gap_penalty (int): The penalty for introducing a gap. Default is -2.

    Returns:
        list: A list of aligned sequences (MSA).
    """
    n = len(sequences)
    distance_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            aligned_a, aligned_b = needleman_wunsch(sequences[i], sequences[j], match_score, mismatch_score, gap_penalty)
            score = sum([match_score if a == b else mismatch_score if a != '-' and b != '-' else gap_penalty
                         for a, b in zip(aligned_a, aligned_b)])
            distance_matrix[i][j] = distance_matrix[j][i] = -score

    center_index = np.argmin(np.sum(distance_matrix, axis=0))
    center_seq = sequences[center_index]
    msa = [center_seq]

    for i in range(n):
        if i == center_index:
            continue
        aligned_center, aligned_seq = needleman_wunsch(center_seq, sequences[i], match_score, mismatch_score, gap_penalty)
        msa = merge_alignments(msa, aligned_center, aligned_seq)
    return msa


def merge_alignments(msa, aligned_center, aligned_seq):
    """
    Merges a newly aligned sequence into the existing MSA using the aligned center sequence
    as reference for insertion of gaps.

    Args:
        msa (list): Current multiple sequence alignment list.
        aligned_center (str): Aligned center sequence from the latest pairwise alignment.
        aligned_seq (str): Newly aligned sequence to add to the MSA.

    Returns:
        list: Updated MSA with the new sequence included.
    """
    new_msa = []

    for seq in msa:
        new_seq = ""
        index = 0
        for i in range(len(aligned_center)):
            if aligned_center[i] == '-':
                new_seq += '-'
            else:
                new_seq += seq[index]
                index += 1
        new_msa.append(new_seq)

    new_msa.append(aligned_seq)

    max_len = max(len(seq) for seq in new_msa)
    new_msa = [seq.ljust(max_len, '-') for seq in new_msa]

    return new_msa


def calculate_statistics(msa):
    """
    Calculates alignment statistics for a multiple sequence alignment.

    Args:
        msa (list): A list of aligned sequences.

    Returns:
        tuple: A tuple containing:
            - matches (int): Number of fully conserved positions (all same base).
            - mismatches (int): Number of positions with substitutions.
            - gaps (int): Total number of gap characters across all sequences.
            - identity (float): Percentage identity across all columns.
    """
    matches = mismatches = gaps = 0
    length = len(msa[0])
    for i in range(length):
        column = [seq[i] for seq in msa]
        if '-' in column:
            gaps += column.count('-')
        else:
            if all(base == column[0] for base in column):
                matches += 1
            else:
                mismatches += 1
    identity = (matches / length) * 100
    return matches, mismatches, gaps, identity


def draw_alignment_grid(c, msa, start_x, start_y, cell_size=12):
    """
    Draws a colored grid representation of the MSA in a PDF using the ReportLab canvas.

    Args:
        c (canvas.Canvas): ReportLab canvas object.
        msa (list): List of aligned sequences.
        start_x (int): X-coordinate to start drawing the grid.
        start_y (int): Y-coordinate to start drawing the grid.
        cell_size (int): Size of each cell in the grid. Default is 12.
    """

    base_colors = {
        'A': colors.lightgreen,
        'C': colors.lightblue,
        'G': colors.khaki,
        'T': colors.salmon,
        '-': colors.whitesmoke
    }

    rows = len(msa)
    cols = len(msa[0])

    for row in range(rows):
        for col in range(cols):
            base = msa[row][col]
            color = base_colors.get(base.upper(), colors.white)
            x = start_x + col * cell_size
            y = start_y - row * cell_size
            c.setFillColor(color)
            c.rect(x, y, cell_size, cell_size, fill=1, stroke=0)
            c.setFillColor(colors.black)
            c.setFont("Helvetica", 8)
            c.drawCentredString(x + cell_size / 2, y + 2, base)


def save_to_pdf(msa, matches, mismatches, gaps, identity, match_score, mismatch_score, gap_penalty):
    """
    Saves the multiple sequence alignment result and summary to a formatted PDF.

    Args:
        msa (list): List of aligned sequences.
        matches (int): Number of matching positions.
        mismatches (int): Number of mismatching positions.
        gaps (int): Total number of gap characters.
        identity (float): Percentage identity between sequences.
        match_score (int): Score for matches.
        mismatch_score (int): Penalty for mismatches.
        gap_penalty (int): Penalty for gaps.

    Returns:
        None
    """
    pdf_filename = "msa_report.pdf"
    c = canvas.Canvas(pdf_filename, pagesize=letter)
    width, height = letter
    margin = 50
    y = height - margin

    c.setFont("Helvetica-Bold", 16)
    c.drawString(100, y, "Center Star Multiple Sequence Alignment Report")
    y -= 50

    c.setFont("Helvetica-Bold", 12)
    c.drawString(100, y, "Alignment Parameters")
    y -= 20
    c.setFont("Helvetica", 12)
    c.drawString(100, y, f"• Match score:        {match_score}")
    y -= 20
    c.drawString(100, y, f"• Mismatch penalty:   {mismatch_score}")
    y -= 20
    c.drawString(100, y, f"• Gap penalty:        {gap_penalty}")
    y -= 30

    c.setFont("Helvetica-Bold", 12)
    c.drawString(100, y, "Alignment Summary")
    y -= 20
    c.setFont("Helvetica", 12)
    c.drawString(100, y, f"• Matches:            {matches}")
    y -= 20
    c.drawString(100, y, f"• Mismatches:         {mismatches}")
    y -= 20
    c.drawString(100, y, f"• Gaps:               {gaps}")
    y -= 20
    c.drawString(100, y, f"• Identity:           {identity:.2f}%")
    y -= 30

    c.setFont("Helvetica-Bold", 12)
    c.drawString(100, y, "Multiple Sequence Alignment (Text Format)")
    y -= 20
    max_line_width = 80
    for i in range(0, len(msa[0]), max_line_width):
        for seq in msa:
            c.drawString(100, y, seq[i:i + max_line_width])
            y -= 20
        y -= 10

    if y < 200:
        c.showPage()
        y = height - margin

    y -= 20
    c.setFont("Helvetica-Bold", 12)
    c.drawString(100, y, "Colored Alignment Grid")
    y -= 20

    draw_alignment_grid(c, msa, start_x=100, start_y=y)
    c.save()
    print(f"PDF saved as {pdf_filename}")


class AlignmentApp:
    """
    GUI application for performing Center-Star Multiple Sequence Alignment (MSA) using Tkinter.

    This class provides a user-friendly interface for inputting DNA sequences (manually or via FASTA file),
    setting alignment parameters (match score, mismatch penalty, gap penalty), and running the MSA process.
    Results are displayed in the GUI and saved to a PDF report.
    """

    def __init__(self, root):
        """
        Initializes the GUI elements and layout for the MSA application.

        Args:
            root (tk.Tk): The root Tkinter window.
        """
        self.root = root
        self.root.title("Multiple Sequence Alignment")
        self.center_window(450, 600)

        # Adding a prompt for the user to input sequences
        tk.Label(root, text="Enter DNA sequences (one per line):").grid(row=0, column=0, padx=10, pady=5, columnspan=2, sticky="w")

        self.seq_text = tk.Text(root, height=10, width=50)
        self.seq_text.grid(row=1, column=0, padx=10, pady=10, columnspan=2)

        # File upload button for FASTA file
        self.upload_button = tk.Button(root, text="Upload FASTA File", command=self.upload_file)
        self.upload_button.grid(row=2, column=0, padx=10, pady=5)

        # Scoring parameters
        self.match_score = tk.IntVar(value=1)
        self.mismatch_penalty = tk.IntVar(value=0)
        self.gap_penalty = tk.IntVar(value=-1)

        tk.Label(root, text="Match Score:").grid(row=3, column=0, padx=10, pady=5, sticky="w")
        tk.Entry(root, textvariable=self.match_score).grid(row=3, column=1, padx=10, pady=5)

        tk.Label(root, text="Mismatch Penalty:").grid(row=4, column=0, padx=10, pady=5, sticky="w")
        tk.Entry(root, textvariable=self.mismatch_penalty).grid(row=4, column=1, padx=10, pady=5)

        tk.Label(root, text="Gap Penalty:").grid(row=5, column=0, padx=10, pady=5, sticky="w")
        tk.Entry(root, textvariable=self.gap_penalty).grid(row=5, column=1, padx=10, pady=5)

        # Result text box
        self.result_text = tk.Text(root, height=10, width=50)
        self.result_text.grid(row=6, column=0, columnspan=2, padx=10, pady=10)

        # Buttons for running alignment and clearing fields
        self.run_button = tk.Button(root, text="Run Alignment", command=self.run_alignment)
        self.run_button.grid(row=7, column=0, padx=10, pady=10)

        self.clear_button = tk.Button(root, text="Clear", command=self.clear_fields)
        self.clear_button.grid(row=7, column=1, padx=10, pady=10)

    def center_window(self, width, height):
        """
        Centers the GUI window on the screen.

        Args:
            width (int): Width of the window.
            height (int): Height of the window.
        """
        # Get the screen width and height
        screen_width = self.root.winfo_screenwidth()
        screen_height = self.root.winfo_screenheight()

        # Calculate the position to center the window
        position_top = int(screen_height / 2 - height / 2)
        position_left = int(screen_width / 2 - width / 2)

        # Set the geometry of the window
        self.root.geometry(f'{width}x{height}+{position_left}+{position_top}')

    def upload_file(self):
        """
        Opens a file dialog for the user to select a FASTA file.
        Loads valid DNA sequences from the file into the input text box.
        Displays an error message if the file is invalid or contains fewer than two sequences.
        """
        # Open a file dialog to select a FASTA file
        file_path = filedialog.askopenfilename(filetypes=[("FASTA files", ("*.fasta", "*.fa")), ("All files", "*.*")])
        if file_path:
            try:
                with open(file_path, 'r') as f:
                    sequences = [str(record.seq) for record in SeqIO.parse(f, "fasta")]
                    if len(sequences) < 2:
                        raise ValueError("FASTA file must contain at least two sequences.")
                    # Automatically populate the text field with sequences
                    self.seq_text.delete(1.0, tk.END)
                    self.seq_text.insert(tk.END, "\n".join(sequences))
            except Exception as e:
                messagebox.showerror("Error", f"An error occurred while loading the file: {e}")

    def run_alignment(self):
        """
        Runs the multiple sequence alignment based on the input sequences and scoring parameters.
        Displays the aligned sequences and statistics in the result box.
        Saves the results in a styled PDF report.
        Shows error messages if the input is invalid or alignment fails.
        """
        seqs = self.seq_text.get("1.0", tk.END).strip().splitlines()
        try:
            sequences = parse_sequences(seqs=seqs)
            match_score = self.match_score.get()
            mismatch_penalty = self.mismatch_penalty.get()
            gap_penalty = self.gap_penalty.get()

            msa = center_star_msa(sequences, match_score, mismatch_penalty, gap_penalty)
            matches, mismatches, gaps, identity = calculate_statistics(msa)

            result_str = f"Matches: {matches}\nMismatches: {mismatches}\nGaps: {gaps}\nIdentity: {identity:.2f}%\n\n"
            result_str += "\n".join(msa)
            self.result_text.delete(1.0, tk.END)
            self.result_text.insert(tk.END, result_str)

            # Save to PDF and show confirmation
            save_to_pdf(msa, matches, mismatches, gaps, identity, match_score, mismatch_penalty, gap_penalty)
            messagebox.showinfo("Success", "PDF file has been created and saved.")
        except Exception as e:
            messagebox.showerror("Error", str(e))

    def clear_fields(self):
        """
        Clears the sequence input and result output text fields.
        """
        self.seq_text.delete(1.0, tk.END)
        self.result_text.delete(1.0, tk.END)


if __name__ == "__main__":
    root = tk.Tk()
    app = AlignmentApp(root)
    root.mainloop()
