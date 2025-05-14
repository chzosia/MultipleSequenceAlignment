# Multiple Sequence Alignment (MSA) with Graphical User Interface (GUI)

This project implements the **Multiple Sequence Alignment (MSA)** using the center star method. The program allows users to load DNA sequences via a user-friendly interface, perform alignment, and view the results, including matrix visualization and alignment summaries in a PDF file.

## Features

- User-friendly GUI: No command-line experience required. Everything is done through an easy-to-use interface.
- Input via FASTA file or direct text input: You can load DNA sequences from a FASTA file or enter them manually.
- Validation of DNA sequences: The program ensures that only valid DNA sequences (A, C, G, T) are entered.
- Multiple Sequence Alignment: Performs multiple sequence alignment.
- Matrix Visualization: Displays a matrix showing the alignment results.
- Auto-generated PDF report: A PDF report is generated containing an alignment summary and graphical representation.
---

## How to Run

### Step 1: Launch the Application
Run the program by launching it from your IDE. The GUI will open up, allowing you to interact with the application.
You can also start it from command line:

```bash
python MSA.py
```
### Step 2: Load DNA Sequences

**From a FASTA File:** Click the "Upload FASTA File" button to select and load a FASTA file containing your DNA sequences.

**Manual Input:** Alternatively, you can manually input your DNA sequences into the provided text fields. Ensure that each sequence is entered in the proper format, with only A, C, G, T characters.

### Step 3: Align Sequences

Once your sequences are loaded, click the **"Run Alignment"** button. The program will perform the multiple sequence alignment and display the results in the GUI.

### Step 4: View Alignment Results

**Textual Alignment:** After alignment, the program will display an easily readable visualization of the aligned sequences. It will also display the match, mismatch, gap and identity statistics.

### Step 5: Check a generated PDF Report

A customized PDF report will be generated and saved on your computer. 
The PDF will include:
  - the alignment parameters
  - summary of alignment statistics 
  - a colored alignment matrix.
