import sys
import requests
from Bio import Entrez
from Bio import SeqIO
import io

from PyQt5.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout,
                             QHBoxLayout, QPushButton, QTextEdit, QLabel,
                             QFileDialog, QMessageBox, QLineEdit, QGroupBox,
                             QListWidget, QSplitter, QFrame, QDialog, QSizePolicy,
                             QProgressBar)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from PyQt5.QtCore import Qt, QDateTime

from PyQt5.QtGui import QTextDocument
from PyQt5.QtPrintSupport import QPrinter
import os
from datetime import datetime


class NCBIDownloadThread(QThread):
    """Wątek do pobierania danych z NCBI w tle"""
    finished = pyqtSignal(object)
    error = pyqtSignal(str)
    progress = pyqtSignal(str)

    def __init__(self, accession_id, email):
        super().__init__()
        self.accession_id = accession_id
        self.email = email

    def run(self):
        try:
            self.progress.emit("Connecting to NCBI...")
            Entrez.email = self.email

            self.progress.emit(f"Downloading sequence: {self.accession_id}")
            # Pobieranie sekwencji w formacie FASTA
            handle = Entrez.efetch(db="nucleotide", id=self.accession_id, rettype="fasta", retmode="text")
            fasta_data = handle.read()
            handle.close()

            self.progress.emit("Processing data...")
            self.finished.emit(fasta_data)

        except Exception as e:
            self.error.emit(f"Error downloading from NCBI: {str(e)}")


class DNAViewerApp(QMainWindow):
    def __init__(self):
        super().__init__()
        self.initUI()

    def initUI(self):
        self.setWindowTitle('DNA Sequence Viewer with NCBI GenBank')
        self.setGeometry(100, 100, 1200, 700)

        # Central widget
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        # Main layout
        main_layout = QHBoxLayout(central_widget)

        # Left panel for controls
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)

        # Load buttons
        load_buttons_layout = QHBoxLayout()
        self.load_file_button = QPushButton('Load from File')
        self.load_file_button.clicked.connect(self.load_dna_data)
        self.load_ncbi_button = QPushButton('Load from NCBI')
        self.load_ncbi_button.clicked.connect(self.load_from_ncbi)
        load_buttons_layout.addWidget(self.load_file_button)
        load_buttons_layout.addWidget(self.load_ncbi_button)
        left_layout.addLayout(load_buttons_layout)

        # NCBI Download section
        ncbi_group = QGroupBox("NCBI GenBank Download")
        ncbi_layout = QVBoxLayout(ncbi_group)

        # Accession ID input
        accession_layout = QHBoxLayout()
        accession_label = QLabel('Accession ID:')
        self.accession_input = QLineEdit()
        self.accession_input.setPlaceholderText('e.g., NM_001301717.1')
        self.accession_input.setText('NM_001301717.1')  # Przykładowy ID
        accession_layout.addWidget(accession_label)
        accession_layout.addWidget(self.accession_input)
        ncbi_layout.addLayout(accession_layout)

        # Email input (required for NCBI)
        email_layout = QHBoxLayout()
        email_label = QLabel('Email:')
        self.email_input = QLineEdit()
        self.email_input.setPlaceholderText('your.email@example.com')
        self.email_input.setText('your.email@example.com')
        email_layout.addWidget(email_label)
        email_layout.addWidget(self.email_input)
        ncbi_layout.addLayout(email_layout)

        # Progress bar for NCBI download
        self.ncbi_progress = QProgressBar()
        self.ncbi_progress.setVisible(False)
        ncbi_layout.addWidget(self.ncbi_progress)

        # Progress label
        self.ncbi_progress_label = QLabel('')
        ncbi_layout.addWidget(self.ncbi_progress_label)

        left_layout.addWidget(ncbi_group)

        # Search group
        search_group = QGroupBox("Sequence Search")
        search_layout = QVBoxLayout(search_group)

        # Search input
        search_input_layout = QHBoxLayout()
        search_label = QLabel('Search pattern:')
        self.search_input = QLineEdit()
        self.search_input.setPlaceholderText('Enter DNA sequence (e.g., ATGC)')
        search_button = QPushButton('Search')
        search_button.clicked.connect(self.search_sequence)
        search_input_layout.addWidget(search_label)
        search_input_layout.addWidget(self.search_input)
        search_input_layout.addWidget(search_button)
        search_layout.addLayout(search_input_layout)

        # Count button and display
        count_layout = QHBoxLayout()
        count_button = QPushButton('Count Occurrences')
        count_button.clicked.connect(self.count_occurrences)
        self.count_result = QLabel('Count: 0')
        self.count_result.setStyleSheet("QLabel { background-color: #f0f0f0; padding: 5px; border: 1px solid #ccc; }")
        count_layout.addWidget(count_button)
        count_layout.addWidget(self.count_result)
        count_layout.addStretch()
        search_layout.addLayout(count_layout)

        # Plot button
        plot_button = QPushButton('Plot Motif Distribution')
        plot_button.clicked.connect(self.plot_motif_distribution)
        search_layout.addWidget(plot_button)

        # Save results button
        self.save_button = QPushButton('Save Results to CSV')
        self.save_button.clicked.connect(self.save_results_to_csv)
        self.save_button.setEnabled(False)  # Initially disabled
        search_layout.addWidget(self.save_button)

        # Save to PDF button
        self.pdf_button = QPushButton('Save Results to PDF')
        self.pdf_button.clicked.connect(self.save_results_to_pdf)
        self.pdf_button.setEnabled(False)  # Initially disabled
        search_layout.addWidget(self.pdf_button)

        # Separator
        separator = QFrame()
        separator.setFrameShape(QFrame.HLine)
        separator.setFrameShadow(QFrame.Sunken)
        search_layout.addWidget(separator)

        # Search results
        results_label = QLabel('Search results:')
        search_layout.addWidget(results_label)

        self.results_list = QListWidget()
        self.results_list.itemClicked.connect(self.highlight_match)
        search_layout.addWidget(self.results_list)

        left_layout.addWidget(search_group)
        left_layout.addStretch()

        # Right panel for data display
        right_panel = QWidget()
        right_layout = QVBoxLayout(right_panel)

        # Sequence ID display
        id_label = QLabel('Sequence ID:')
        right_layout.addWidget(id_label)

        self.id_display = QTextEdit()
        self.id_display.setMaximumHeight(50)
        self.id_display.setReadOnly(True)
        right_layout.addWidget(self.id_display)

        # Sequence display
        sequence_label = QLabel('DNA Sequence:')
        right_layout.addWidget(sequence_label)

        self.sequence_display = QTextEdit()
        self.sequence_display.setReadOnly(True)
        self.sequence_display.setMaximumHeight(100)  # Zmniejszona wysokość
        right_layout.addWidget(self.sequence_display)

        # Splitter to allow resizing between sequence and plot
        splitter = QSplitter(Qt.Vertical)

        # Create a container for the plot with proper size policy
        plot_container = QWidget()
        plot_layout = QVBoxLayout(plot_container)
        plot_layout.setContentsMargins(0, 0, 0, 0)

        # Plot label
        plot_label = QLabel('Motif Distribution:')
        plot_layout.addWidget(plot_label)

        # Create matplotlib figure and canvas
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.canvas.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.toolbar = NavigationToolbar(self.canvas, self)

        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)

        # Add widgets to splitter
        splitter.addWidget(self.sequence_display)
        splitter.addWidget(plot_container)

        # Set initial sizes (sequence display takes 30%, plot takes 70%)
        splitter.setSizes([100, 500])

        right_layout.addWidget(splitter)

        # Add panels to main layout
        main_layout.addWidget(left_panel, 1)
        main_layout.addWidget(right_panel, 3)

    def load_dna_data(self):
        # Open file dialog to select DNA sequence file
        file_path, _ = QFileDialog.getOpenFileName(
            self, 'Open DNA Sequence File', '',
            'Text Files (*.txt);;FASTA Files (*.fasta *.fa);;All Files (*)'
        )

        if file_path:
            try:
                with open(file_path, 'r') as file:
                    content = file.read().strip()

                self.parse_sequence_content(content, f"File: {os.path.basename(file_path)}")

            except Exception as e:
                QMessageBox.critical(self, 'Error', f'Could not read file: {str(e)}')

    def load_from_ncbi(self):
        """Load sequence from NCBI GenBank using accession ID"""
        accession_id = self.accession_input.text().strip()
        email = self.email_input.text().strip()

        if not accession_id:
            QMessageBox.warning(self, 'Error', 'Please enter an accession ID')
            return

        if not email or '@' not in email:
            QMessageBox.warning(self, 'Error', 'Please enter a valid email address')
            return

        # Disable button during download
        self.load_ncbi_button.setEnabled(False)
        self.ncbi_progress.setVisible(True)
        self.ncbi_progress_label.setText('Starting download...')

        # Create and start download thread
        self.ncbi_thread = NCBIDownloadThread(accession_id, email)
        self.ncbi_thread.finished.connect(self.on_ncbi_download_finished)
        self.ncbi_thread.error.connect(self.on_ncbi_download_error)
        self.ncbi_thread.progress.connect(self.on_ncbi_progress)
        self.ncbi_thread.start()

    def on_ncbi_progress(self, message):
        """Update progress during NCBI download"""
        self.ncbi_progress_label.setText(message)

    def on_ncbi_download_finished(self, fasta_data):
        """Handle successful NCBI download"""
        try:
            self.parse_sequence_content(fasta_data, f"NCBI: {self.accession_input.text()}")
            self.ncbi_progress_label.setText('Download completed successfully!')
            QMessageBox.information(self, 'Success', 'Sequence downloaded successfully from NCBI GenBank')
        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Error processing NCBI data: {str(e)}')
        finally:
            self.load_ncbi_button.setEnabled(True)
            self.ncbi_progress.setVisible(False)

    def on_ncbi_download_error(self, error_message):
        """Handle NCBI download error"""
        QMessageBox.critical(self, 'Error', error_message)
        self.load_ncbi_button.setEnabled(True)
        self.ncbi_progress.setVisible(False)
        self.ncbi_progress_label.setText('')

    def parse_sequence_content(self, content, source_info):
        """Parse sequence content from either file or NCBI data"""
        lines = content.split('\n')

        if len(lines) > 0:
            # First line is typically the ID line (starts with '>')
            if lines[0].startswith('>'):
                sequence_id = lines[0][1:].strip()  # Remove '>' and trim
                sequence = ''.join(lines[1:]).replace(' ', '').upper()
                # Add source info to ID
                sequence_id = f"{sequence_id} | {source_info}"
            else:
                # If no '>', treat first line as ID and rest as sequence
                sequence_id = lines[0].strip()
                sequence = ''.join(lines[1:]).replace(' ', '').upper()
                sequence_id = f"{sequence_id} | {source_info}"

            # Validate DNA sequence
            if not self.validate_dna_sequence(sequence):
                QMessageBox.warning(self, 'Warning',
                                    'Sequence contains non-DNA characters. Only A, T, C, G will be processed.')

            # Display the data
            self.id_display.setText(sequence_id)
            self.sequence_display.setText(sequence)
            self.current_sequence = sequence
            self.clear_search_results()

            # Clear the plot when new data is loaded
            self.figure.clear()
            self.canvas.draw()

        else:
            QMessageBox.warning(self, 'Error', 'No data found')

    def validate_dna_sequence(self, sequence):
        """Validate if sequence contains only DNA characters"""
        valid_bases = {'A', 'T', 'C', 'G', 'N', ' '}
        return all(base.upper() in valid_bases for base in sequence)

    def search_sequence(self):
        if not hasattr(self, 'current_sequence'):
            QMessageBox.warning(self, 'Error', 'Please load a DNA sequence first')
            return

        search_pattern = self.search_input.text().strip().upper()

        if not search_pattern:
            QMessageBox.warning(self, 'Error', 'Please enter a search pattern')
            return

        # Validate DNA sequence (only A, T, C, G characters)
        valid_bases = {'A', 'T', 'C', 'G'}
        if not all(base in valid_bases for base in search_pattern):
            QMessageBox.warning(self, 'Error', 'Invalid DNA sequence. Only A, T, C, G characters are allowed.')
            return

        # Perform search
        sequence = self.current_sequence
        matches = []

        # Find all occurrences
        start = 0
        while True:
            pos = sequence.find(search_pattern, start)
            if pos == -1:
                break
            matches.append((pos, pos + len(search_pattern)))
            start = pos + 1

        # Display results
        self.results_list.clear()

        if not matches:
            self.results_list.addItem('No matches found')
            self.save_button.setEnabled(False)  # Disable save if no matches
            return

        for i, (start, end) in enumerate(matches):
            self.results_list.addItem(f'Match {i + 1}: Position {start}-{end - 1}')

        self.matches = matches
        self.search_pattern = search_pattern
        self.save_button.setEnabled(True)  # Enable save button
        self.pdf_button.setEnabled(True)  # Enable PDF button

    def highlight_match(self, item):
        if not hasattr(self, 'matches') or not self.matches:
            return

        # Get the selected match index
        try:
            index = self.results_list.row(item)
            if index >= len(self.matches):
                return

            start, end = self.matches[index]

            # Highlight the match in the sequence display
            cursor = self.sequence_display.textCursor()

            # Set selection
            cursor.setPosition(0)
            cursor.movePosition(cursor.Right, cursor.MoveAnchor, start)
            cursor.movePosition(cursor.Right, cursor.KeepAnchor, end - start)

            self.sequence_display.setTextCursor(cursor)
            self.sequence_display.setFocus()

            # Apply highlighting (optional - you can customize this)
            extra_selections = []
            selection = QTextEdit.ExtraSelection()
            selection.cursor = cursor
            selection.format.setBackground(Qt.yellow)
            selection.format.setProperty(QTextEdit.FullWidthSelection, True)
            extra_selections.append(selection)

            self.sequence_display.setExtraSelections(extra_selections)

        except Exception as e:
            print(f"Error highlighting match: {e}")

    def count_occurrences(self):
        if not hasattr(self, 'current_sequence'):
            QMessageBox.warning(self, 'Error', 'Please load a DNA sequence first')
            return

        search_pattern = self.search_input.text().strip().upper()

        if not search_pattern:
            QMessageBox.warning(self, 'Error', 'Please enter a search pattern')
            return

        # Validate DNA sequence (only A, T, C, G characters)
        valid_bases = {'A', 'T', 'C', 'G'}
        if not all(base in valid_bases for base in search_pattern):
            QMessageBox.warning(self, 'Error', 'Invalid DNA sequence. Only A, T, C, G characters are allowed.')
            return

        # Count occurrences
        sequence = self.current_sequence
        count = sequence.count(search_pattern)

        # Update count display
        self.count_result.setText(f'Count: {count}')

        # Also perform search to update results list
        self.search_sequence()

        # Show message with count result
        QMessageBox.information(self, 'Count Result',
                                f'The pattern "{search_pattern}" appears {count} times in the sequence.')

    def plot_motif_distribution(self):
        if not hasattr(self, 'current_sequence'):
            QMessageBox.warning(self, 'Error', 'Please load a DNA sequence first')
            return

        if not hasattr(self, 'matches') or not self.matches:
            QMessageBox.warning(self, 'Error', 'Please search for a pattern first')
            return

            # Clear previous plot
        self.figure.clear()

        # Create subplots
        ax1 = self.figure.add_subplot(211)  # Top plot: distribution histogram
        ax2 = self.figure.add_subplot(212)  # Bottom plot: cumulative distribution

        # Get sequence length and positions
        seq_len = len(self.current_sequence)
        positions = [start for start, end in self.matches]

        # Plot 1: Histogram of motif distribution
        n_bins = min(20, max(5, seq_len // 100))  # Adaptive bin count
        ax1.hist(positions, bins=n_bins, color='skyblue', edgecolor='black', alpha=0.7)
        ax1.set_title(f'Distribution of "{self.search_pattern}" motifs in DNA sequence')
        ax1.set_xlabel('Position in sequence')
        ax1.set_ylabel('Frequency')
        ax1.grid(True, alpha=0.3)

        # Add vertical lines for each motif position
        for pos in positions:
            ax1.axvline(x=pos, color='red', linestyle='--', alpha=0.3, linewidth=0.5)

        # Plot 2: Cumulative distribution
        if len(positions) > 0:
            sorted_positions = sorted(positions)
            cumulative_counts = range(1, len(sorted_positions) + 1)
            ax2.plot(sorted_positions, cumulative_counts, 'b-', linewidth=2, marker='o', markersize=3)
            ax2.fill_between(sorted_positions, cumulative_counts, alpha=0.3, color='blue')
        ax2.set_xlabel('Position in sequence')
        ax2.set_ylabel('Cumulative Count')
        ax2.set_title('Cumulative Distribution of Motifs')
        ax2.grid(True, alpha=0.3)

        # Add statistics text
        stats_text = (
            f'Total motifs: {len(self.matches)}\n'
            f'Sequence length: {seq_len} bp\n'
            f'Density: {len(self.matches) / seq_len * 1000:.2f} motifs/kb\n'
            f'Mean spacing: {seq_len / len(self.matches):.1f} bp' if len(self.matches) > 0 else 'N/A'
        )

        ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                 fontsize=9)

        # Adjust layout to prevent clipping
        self.figure.tight_layout()

        # Refresh canvas
        self.canvas.draw()

    def clear_search_results(self):
        self.results_list.clear()
        if hasattr(self, 'matches'):
            del self.matches
        if hasattr(self, 'search_pattern'):
            del self.search_pattern
        # Clear any highlighting
        self.sequence_display.setExtraSelections([])
        # Disable save button
        self.save_button.setEnabled(False)
        self.pdf_button.setEnabled(False)

    def save_results_to_csv(self):
        if not hasattr(self, 'current_sequence'):
            QMessageBox.warning(self, 'Error', 'Please load a DNA sequence first')
            return

        if not hasattr(self, 'matches') or not self.matches:
            QMessageBox.warning(self, 'Error', 'Please search for a pattern first')
            return

        # Open file dialog to choose save location
        file_path, _ = QFileDialog.getSaveFileName(
            self, 'Save Results to CSV', '',
            'CSV Files (*.csv);;All Files (*)'
        )

        if not file_path:
            return  # User cancelled

        try:
            with open(file_path, 'w', newline='', encoding='utf-8') as csvfile:
                import csv
                writer = csv.writer(csvfile)

                # Write header
                writer.writerow(['DNA Sequence Analysis Results'])
                writer.writerow([])

                # Write sequence info
                writer.writerow(['Sequence ID:', self.id_display.toPlainText()])
                writer.writerow(['Sequence Length:', len(self.current_sequence)])
                writer.writerow(['Search Pattern:', self.search_pattern])
                writer.writerow(['Total Matches:', len(self.matches)])
                writer.writerow([])

                # Write matches table header
                writer.writerow(['Match #', 'Start Position', 'End Position', 'Sequence Segment'])

                # Write each match with sequence context
                for i, (start, end) in enumerate(self.matches):
                    # Get the matched sequence segment
                    sequence_segment = self.current_sequence[start:end]
                    writer.writerow([i + 1, start, end - 1, sequence_segment])

                writer.writerow([])

                # Write distribution data for plotting
                writer.writerow(['Position Distribution Data'])
                writer.writerow(['Position', 'Motif Presence'])
                positions = [start for start, end in self.matches]
                for pos in positions:
                    writer.writerow([pos, 1])

                writer.writerow([])
                writer.writerow(['Generated on:', QDateTime.currentDateTime().toString('yyyy-MM-dd hh:mm:ss')])

            QMessageBox.information(self, 'Success', f'Results saved to {file_path}')

        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Could not save file: {str(e)}')

    def save_results_to_pdf(self):
        if not hasattr(self, 'current_sequence'):
            QMessageBox.warning(self, 'Error', 'Please load a DNA sequence first')
            return

        if not hasattr(self, 'matches') or not self.matches:
            QMessageBox.warning(self, 'Error', 'Please search for a pattern first')
            return

        # Open file dialog to choose save location
        file_path, _ = QFileDialog.getSaveFileName(
            self, 'Save Results to PDF', '',
            'PDF Files (*.pdf);;All Files (*)'
        )

        if not file_path:
            return  # User cancelled

        try:
            # Create printer and set properties
            printer = QPrinter(QPrinter.HighResolution)
            printer.setOutputFormat(QPrinter.PdfFormat)
            printer.setOutputFileName(file_path)
            printer.setPageSize(QPrinter.A4)
            printer.setFullPage(True)

            # Create document
            document = QTextDocument()

            # Prepare HTML content
            html_content = self.generate_pdf_content()
            document.setHtml(html_content)

            # Print to PDF
            document.print_(printer)

            QMessageBox.information(self, 'Success', f'Results saved to PDF: {file_path}')

        except Exception as e:
            QMessageBox.critical(self, 'Error', f'Could not save PDF: {str(e)}')

    def generate_pdf_content(self):
        """Generate HTML content for PDF report"""

        # Save current plot to temporary image
        plot_filename = "temp_plot.png"
        self.figure.savefig(plot_filename, dpi=300, bbox_inches='tight')

        # Get current date and time
        current_time = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Generate matches table rows
        matches_rows = ""
        for i, (start, end) in enumerate(self.matches):
            sequence_segment = self.current_sequence[start:end]
            matches_rows += f"""
            <tr>
                <td>{i + 1}</td>
                <td>{start}</td>
                <td>{end - 1}</td>
                <td>{sequence_segment}</td>
            </tr>
            """

        # Generate additional statistics for the report
        seq_len = len(self.current_sequence)
        motif_density = len(self.matches) / seq_len * 1000 if seq_len > 0 else 0
        sequence_coverage = (sum(len(self.search_pattern) for _ in self.matches) / seq_len * 100) if seq_len > 0 else 0
        mean_spacing = seq_len / len(self.matches) if len(self.matches) > 0 else 0

        html = f"""
        <html>
        <head>
            <style>
                body {{ font-family: Arial, sans-serif; margin: 20px; }}
                h1 {{ color: #2c3e50; border-bottom: 2px solid #3498db; }}
                h2 {{ color: #34495e; }}
                .header {{ background-color: #ecf0f1; padding: 15px; border-radius: 5px; }}
                .section {{ margin: 20px 0; }}
                table {{ width: 100%; border-collapse: collapse; margin: 10px 0; }}
                th, td {{ border: 1px solid #bdc3c7; padding: 8px; text-align: left; }}
                th {{ background-color: #3498db; color: white; }}
                tr:nth-child(even) {{ background-color: #f2f2f2; }}
                .stats {{ background-color: #e8f4f8; padding: 10px; border-radius: 5px; }}
                .footer {{ margin-top: 30px; font-size: 12px; color: #7f8c8d; }}
                .plot-container {{ text-align: center; margin: 20px 0; }}
                .plot-image {{ max-width: 100%; height: auto; border: 1px solid #ddd; }}
            </style>
        </head>
        <body>
            <div class="header">
                <h1>DNA Sequence Analysis Report</h1>
                <p>Generated on: {current_time}</p>
            </div>

            <div class="section">
                <h2>Sequence Information</h2>
                <p><strong>Sequence ID:</strong> {self.id_display.toPlainText()}</p>
                <p><strong>Sequence Length:</strong> {seq_len} bp</p>
                <p><strong>Search Pattern:</strong> {self.search_pattern}</p>
                <p><strong>Pattern Length:</strong> {len(self.search_pattern)} bp</p>
            </div>

            <div class="section stats">
                <h2>Analysis Statistics</h2>
                <p><strong>Total Matches Found:</strong> {len(self.matches)}</p>
                <p><strong>Motif Density:</strong> {motif_density:.2f} motifs/kb</p>
                <p><strong>Sequence Coverage:</strong> {sequence_coverage:.2f}%</p>
                <p><strong>Mean Spacing:</strong> {mean_spacing:.1f} bp</p>
                <p><strong>GC Content of Pattern:</strong> {(self.search_pattern.count('G') + self.search_pattern.count('C')) / len(self.search_pattern) * 100:.1f}%</p>
            </div>

            <div class="section">
                <h2>Motif Distribution Analysis</h2>
                <div class="plot-container">
                    <img src="{plot_filename}" class="plot-image" alt="Motif Distribution Plot" />
                    <p><em>Figure 1: Distribution of '{self.search_pattern}' motifs in the DNA sequence. Top: Frequency distribution histogram. Bottom: Cumulative distribution of motif positions.</em></p>
                </div>
            </div>

            <div class="section">
                <h2>Detailed Match Results</h2>
                <table>
                    <tr>
                        <th>Match #</th>
                        <th>Start Position</th>
                        <th>End Position</th>
                        <th>Sequence Segment</th>
                    </tr>
                    {matches_rows}
                </table>
            </div>

            <div class="section">
                <h2>Position Data</h2>
                <p><strong>Motif Positions:</strong> {', '.join(str(start) for start, end in self.matches[:10])}{'...' if len(self.matches) > 10 else ''}</p>
                <p><strong>Total positions analyzed:</strong> {seq_len} bp</p>
            </div>

            <div class="footer">
                <p>Generated by DNA Sequence Viewer Application</p>
                <p>Total matches: {len(self.matches)} | Sequence length: {seq_len} bp | Pattern: {self.search_pattern}</p>
            </div>
        </body>
        </html>
        """

        # Clean up temporary file
        try:
            os.remove(plot_filename)
        except:
            pass

        return html


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = DNAViewerApp()
    window.show()
    sys.exit(app.exec_())