# RUN DIMPLE
# script for GUI

import ast
import logging
import tkinter as tk
import warnings
from datetime import datetime
from io import StringIO
from tkinter import filedialog, messagebox

from DIMPLE.DIMPLE import (
    addgene,
    align_genevariation,
    generate_DMS_fragments,
    post_qc,
    print_all,
    switch_fragmentsize,
)
from DIMPLE.run_settings import (
    DEFAULT_GUI_RANDOM_SEED,
    DimpleRuntimeConfig,
    apply_barcode_start,
    apply_handle,
    apply_instance_settings,
    apply_random_seed,
    apply_restriction_settings,
    apply_runtime_policies,
    compute_overlaps_and_maxfrag,
    configure_dimple_logging,
    normalize_avoid_list,
    resolve_codon_usage,
    validate_insertions,
)
from DIMPLE.utilities import parse_custom_mutations

log_file = "logs/DIMPLE-{:%Y-%m-%d-%s}.log".format(datetime.now())
configure_dimple_logging(log_file)
logger = logging.getLogger(__name__)
logger.info("Started logging")

AMINO_ACIDS = [
    "Cys",
    "Asp",
    "Ser",
    "Gln",
    "Met",
    "Asn",
    "Pro",
    "Lys",
    "Thr",
    "Phe",
    "Ala",
    "Gly",
    "Ile",
    "Leu",
    "His",
    "Arg",
    "Trp",
    "Val",
    "Glu",
    "Tyr",
]


def run():
    # A fresh DimpleRuntimeConfig is constructed every invocation, which gives
    # this run its own barcode pool (config.barcode_f). Barcodes therefore reset
    # between runs without any explicit teardown -- see upstream issue #21.
    runtime_config = DimpleRuntimeConfig()
    # Check that a mutation type is selected
    if not any(
        [
            app.delete.get(),
            app.insert.get(),
            app.include_substitutions.get(),
            app.dis.get(),
        ]
    ):
        app.output_text.insert(tk.END, "Error: You must select a mutation type.\n")
        messagebox.showerror("Python Error", "Error: You must select a mutation type.")
        raise ValueError("You must select a mutation type")

    # Check if gene file is set
    if app.geneFile is None:
        app.output_text.insert(tk.END, "Error: No gene file selected.\n")
        messagebox.showerror("Python Error", "Error: No gene file selected.")
        raise ValueError("No gene file selected for input.")

    # Set working dir from gene file if not set
    if app.wDir is None:
        app.wDir = app.geneFile.rsplit("/", 1)[0] + "/"
        app.output_text.insert(tk.END, f"Working directory set to {app.wDir}\n")

    try:
        apply_handle(app.handle.get(), config=runtime_config)
    except ValueError as exc:
        app.output_text.insert(tk.END, "Error: Genetic handle contains non-nucleic acid bases\n")
        raise ValueError(str(exc)) from exc

    deletions_for_overlap = False
    if app.delete.get():
        deletions_for_overlap = [int(x) for x in app.deletions.get().split(",")]

    frag_raw = app.fragmentLen.get().strip()
    fragment_len = 0 if frag_raw == "" or frag_raw.lower() == "auto" else int(frag_raw)

    overlap_l, overlap_r = compute_overlaps_and_maxfrag(
        int(app.oligoLen.get()),
        fragment_len,
        int(app.overlap.get()),
        deletions_for_overlap,
        logger=logger,
        config=runtime_config,
    )

    resolve_codon_usage(app.codon_usage, config=runtime_config)

    apply_barcode_start(int(app.barcode_start.get()), config=runtime_config)

    try:
        apply_restriction_settings(app.restriction_sequence.get(), config=runtime_config)
    except (ValueError, IndexError) as exc:
        message = (
            f"Restriction sequence {app.restriction_sequence.get()!r} is not valid. "
            "Expected format like 'CGTCTC(G)1/5' or a named enzyme ('BsaI', 'BsmBI'). "
            f"({exc})"
        )
        app.output_text.insert(tk.END, f"Error: {message}\n")
        messagebox.showerror("Invalid Restriction Sequence", message)
        return

    if runtime_config.enzyme is not None:
        app.output_text.insert(tk.END, f"Using restriction enzyme {runtime_config.enzyme}\n")
    app.output_text.insert(tk.END, f"Using restriction sequence {runtime_config.cutsite}\n")

    normalize_avoid_list(app.avoid_sequence.get(), logger=logger, config=runtime_config)
    apply_runtime_policies(
        dms=bool(app.include_substitutions.get()),
        stop_codon=bool(app.stop.get()),
        make_double=bool(app.make_double.get()),
        maximize_nucleotide_change=bool(app.max_mutations.get()),
        non_interactive=False,
        preferred_orf_index=None,
        link_policy="prompt",
        breaksite_change_policy="prompt",
        config=runtime_config,
    )

    if app.delete.get() == 0:
        deletions = False
    else:
        deletions = [int(x) for x in app.deletions.get().split(",")]
        if any(x % 3 != 0 for x in deletions):
            message = (
                "Warning: chosen deletions will not maintain reading frame. Remember "
                "that deletions are in nucleotides, not codons."
            )
            app.output_text.insert(tk.END, message + "\n")
            warnings.warn(message)
    if app.insert.get() == 0:
        insertions = False
    else:
        insertions = app.insertions.get().split(",")
        if any(len(x) % 3 != 0 for x in insertions):
            message = (
                "Warning: chosen insertions will not maintain reading frame. Remember "
                "that insertions are in nucleotides, not codons."
            )
            app.output_text.insert(tk.END, message + "\n")
            warnings.warn(message)
        try:
            validate_insertions(insertions, runtime_config)
        except ValueError as exc:
            app.output_text.insert(tk.END, f"Error: {exc}\n")
            raise

    apply_random_seed(DEFAULT_GUI_RANDOM_SEED, config=runtime_config)

    pool = addgene(app.geneFile, runtime_config)

    apply_instance_settings(
        pool,
        aminoacids=app.substitutions.get().split(","),
        doublefrag=app.doublefrag.get(),
        gene_primer_tm=(
            int(app.melting_temp_low.get()),
            int(app.melting_temp_high.get()),
        ),
        config=runtime_config,
    )
    if app.avoid_breaksites.get():
        pool[0].problemsites = set(int(x) for x in app.custom_mutations.keys())
        # add extras
        if app.avoid_others_list.get() != "":
            pool[0].problemsites.update([int(x) for x in app.avoid_others_list.get().split(",")])
        for i in range(len(pool[0].breaksites)):
            switch_fragmentsize(pool[0], 1, pool)

    if app.matchSequences.get() == "match":
        align_genevariation(pool)

    logger.info("Generating DMS fragments")
    app.output_text.insert(tk.END, "Generating DMS fragments\n")

    # Generate DMS fragments
    generate_DMS_fragments(
        pool,
        overlap_l,
        overlap_r,
        app.synonymous.get(),
        app.custom_mutations,
        app.include_substitutions.get(),
        insertions,
        deletions,
        app.dis.get(),
        app.wDir,
    )

    # Post QC checks and saving
    app.output_text.insert(tk.END, "Post QC checks and saving\n")
    post_qc(pool)
    print_all(pool, app.wDir)

    logger.info("Finished")
    app.output_text.insert(tk.END, "==========================================\n")
    app.output_text.insert(tk.END, "DONE: DIMPLE run finished successfully.\n")
    app.output_text.insert(tk.END, f"Output written to: {app.wDir}\n")
    app.output_text.insert(tk.END, f"Log file saved to {log_file}\n")
    app.output_text.insert(tk.END, "==========================================\n")
    app.output_text.see(tk.END)
    messagebox.showinfo(
        "DIMPLE",
        f"DIMPLE run finished successfully.\n\nOutput written to:\n{app.wDir}",
    )

    # Output all parameters to log


class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.pack()
        self.create_widgets()

    def create_widgets(self):
        self.winfo_toplevel().title("DIMPLE Deep Indel Missense Programmable Library Engineering")
        self.matchSequences = tk.IntVar()
        self.mutationType = tk.IntVar()
        self.usage = tk.IntVar()
        self.include_substitutions = tk.IntVar()
        self.delete = tk.IntVar()
        self.insert = tk.IntVar()
        self.stop = tk.IntVar()
        self.synonymous = tk.IntVar()
        self.dis = tk.IntVar()
        self.make_double = tk.IntVar()
        self.custom_mutations = {}
        self.doublefrag = tk.IntVar()
        self.avoid_breaksites = tk.IntVar()
        self.max_mutations = tk.IntVar()

        self.wDir_file = tk.Button(self, text="Working Directory", command=self.browse_wDir)
        self.wDir_file.pack()
        self.wDir = None

        self.gene_file = tk.Button(self, text="Target Gene File", command=self.browse_gene)
        self.gene_file.pack()
        self.geneFile = None

        tk.Label(self, text="Oligo Length").pack()
        self.oligoLen = tk.Entry(self, textvariable=tk.StringVar(self, "250"))
        self.oligoLen.pack()

        tk.Label(self, text="Fragment Length ('auto' or integer)").pack()
        self.fragmentLen = tk.Entry(self, textvariable=tk.StringVar(self, "auto"))
        self.fragmentLen.pack()

        tk.Label(self, text="Overlap (bases)").pack()
        self.overlap = tk.Entry(self, textvariable=tk.StringVar(self, "4"))
        self.overlap.pack()

        tk.Label(self, text="Barcode Start position (3000 total available)").pack()
        self.barcode_start = tk.Entry(self, textvariable=tk.StringVar(self, "0"))
        self.barcode_start.pack()

        self.melting_temp_low = tk.Entry(self, textvariable=tk.StringVar(self, "58"))
        self.melting_temp_high = tk.Entry(self, textvariable=tk.StringVar(self, "62"))

        tk.Label(self, text="Type IIS restriction sequence (Do not use N)").pack()
        self.restriction_sequence = tk.Entry(self, textvariable=tk.StringVar(self, "CGTCTC(G)1/5"))
        self.restriction_sequence.pack()

        tk.Label(self, text="Sequences to avoid").pack()
        self.avoid_sequence = tk.Entry(
            self, width=50, textvariable=tk.StringVar(self, "CGTCTC, GGTCTC")
        )
        self.avoid_sequence.pack()

        self.stop_codon = tk.Checkbutton(self, text="Include Stop Codons", variable=self.stop)
        self.stop_codon.pack()

        self.synonymous_check = tk.Checkbutton(
            self, text="Include Synonymous Mutations", variable=self.synonymous
        )
        self.synonymous_check.pack()

        self.doublefrag_check = tk.Checkbutton(
            self, text="Double Fragments per Oligo", variable=self.doublefrag
        )
        self.doublefrag_check.pack()

        self.codon_usage = "human"

        self.custom_codons = {
            "TTT": 0.45,
            "TTC": 0.55,
            "TTA": 0.07,
            "TTG": 0.13,
            "TAT": 0.43,
            "TAC": 0.57,
            "TAA": 0.28,
            "TAG": 0.2,
            "CTT": 0.13,
            "CTC": 0.2,
            "CTA": 0.07,
            "CTG": 0.41,
            "CAT": 0.41,
            "CAC": 0.59,
            "CAA": 0.25,
            "CAG": 0.75,
            "ATT": 0.36,
            "ATC": 0.48,
            "ATA": 0.16,
            "ATG": 1,
            "AAT": 0.46,
            "AAC": 0.54,
            "AAA": 0.42,
            "AAG": 0.58,
            "GTT": 0.18,
            "GTC": 0.24,
            "GTA": 0.11,
            "GTG": 0.47,
            "GAT": 0.46,
            "GAC": 0.54,
            "GAA": 0.42,
            "GAG": 0.58,
            "TCT": 0.18,
            "TCC": 0.22,
            "TCA": 0.15,
            "TCG": 0.06,
            "TGT": 0.45,
            "TGC": 0.55,
            "TGA": 0.52,
            "TGG": 1,
            "CCT": 0.28,
            "CCC": 0.33,
            "CCA": 0.27,
            "CCG": 0.11,
            "CGT": 0.08,
            "CGC": 0.19,
            "CGA": 0.11,
            "CGG": 0.21,
            "ACT": 0.24,
            "ACC": 0.36,
            "ACA": 0.28,
            "ACG": 0.12,
            "AGT": 0.15,
            "AGC": 0.24,
            "AGA": 0.2,
            "AGG": 0.2,
            "GCT": 0.26,
            "GCC": 0.4,
            "GCA": 0.23,
            "GCG": 0.11,
            "GGT": 0.16,
            "GGC": 0.34,
            "GGA": 0.25,
            "GGG": 0.25,
        }

        def openNewWindow():
            newWindow = tk.Toplevel(self)
            newWindow.title("Custom Codon Usage")
            newWindow.geometry("900x400")
            custom_codon = tk.Text(newWindow, width=750, height=400, wrap=tk.WORD)
            custom_codon.pack()
            # if self.codon_usage is not a dictionary then set it to self.custom_codons
            if type(self.codon_usage) is not dict:
                self.codon_usage = self.custom_codons
            s = StringIO()
            print(self.codon_usage, file=s)
            custom_codon.insert(tk.END, s.getvalue())

            def on_closing():
                # set codon usage
                self.codon_usage = ast.literal_eval(custom_codon.get("1.0", "end-1c").strip())
                newWindow.destroy()

            newWindow.protocol("WM_DELETE_WINDOW", on_closing)

        def ecoli_ON():
            self.codon_usage = "ecoli"

        def human_ON():
            self.codon_usage = "human"

        tk.Label(self, text="Codon Usage", font="helvetica 12 underline").pack(pady=5)
        self.ecoli_check = tk.Radiobutton(
            self, text="E. coli", variable=self.usage, value=1, command=ecoli_ON
        )
        self.ecoli_check.pack()
        self.human_check = tk.Radiobutton(
            self, text="Human", variable=self.usage, value=0, command=human_ON
        )
        self.human_check.pack()
        self.custom_codon_button = tk.Button(self, text="Custom Codon Usage", command=openNewWindow)
        self.custom_codon_button.pack()

        tk.Label(self, text="Select Mutations", font="helvetica 12 underline").pack(pady=5)

        self.include_dis = tk.Checkbutton(self, text="Domain Insertion Scan", variable=self.dis)
        self.include_dis.pack()
        self.handle = tk.Entry(
            self, width=50, textvariable=tk.StringVar(self, "AGCGGGAGACCGGGGTCTCTGAGC")
        )
        self.handle.pack()

        self.delete_check = tk.Checkbutton(self, text="List of Deletions", variable=self.delete)
        self.delete_check.pack()
        self.delete_check.deselect()
        self.deletions = tk.Entry(self, width=50, textvariable=tk.StringVar(self, "3,6"))
        self.deletions.pack()

        self.insert_check = tk.Checkbutton(self, text="List of Insertions", variable=self.insert)
        self.insert_check.pack()
        self.insert_check.deselect()
        self.insertions = tk.Entry(
            self, width=50, textvariable=tk.StringVar(self, "GGC,GGCTCT,GGCTCTGGA")
        )
        self.insertions.pack()

        self.include_sub_check = tk.Checkbutton(
            self, text="Deep Mutational Scan", variable=self.include_substitutions
        )
        self.include_sub_check.pack()
        self.include_sub_check.deselect()
        self.max_mut = tk.Checkbutton(
            self,
            text="Maximize Nucleotide Change (2 or more)",
            variable=self.max_mutations,
        )
        self.max_mut.pack()
        self.max_mut.deselect()

        self.substitutions = tk.Entry(
            self,
            width=80,
            textvariable=tk.StringVar(
                self,
                "Cys,Asp,Ser,Gln,Met,Asn,Pro,Lys,Thr,Phe,Ala,Gly,Ile,Leu,His,Arg,Trp,Val,Glu,Tyr",
            ),
        )
        self.substitutions.pack()

        self.double_it = tk.Checkbutton(
            self,
            text="Make Double Mutations (Warning: Limit number of mutations)",
            variable=self.make_double,
        )
        self.double_it.pack()
        self.double_it.deselect()

        self.custom_mutations_button = tk.Button(
            self, text="Custom Mutations", command=self.custom_mutations_input
        )
        self.custom_mutations_button.pack(pady=10)

        self.avoid_custom = tk.Checkbutton(
            self,
            text="Avoid breaksites in Custom Mutations",
            variable=self.avoid_breaksites,
        )
        self.avoid_custom.pack()

        self.avoid_others_list = ""

        # self.matchSequences_check = tk.Checkbutton(
        #     self, text='Match Sequences', variable=self.matchSequences)
        # self.matchSequences_check.pack()

        self.run = tk.Button(self, text="Run DIMPLE", command=run).pack(pady=10)

        # Set up output box
        # Create a frame to hold the Text widget and scrollbar
        output_frame = tk.Frame(self)
        output_frame.pack(fill="both", expand=True, padx=10, pady=10)

        # Create the Text widget
        self.output_text = tk.Text(output_frame, wrap="word", height=20)
        self.output_text.pack(side="left", fill="both", expand=True)

        # Create the vertical scrollbar
        scrollbar = tk.Scrollbar(output_frame, orient="vertical", command=self.output_text.yview)
        scrollbar.pack(side="right", fill="y")

        # Link the scrollbar to the Text widget
        self.output_text.configure(yscrollcommand=scrollbar.set)

    def browse_wDir(self):
        self.wDir = filedialog.askdirectory(title="Select a File")
        if self.wDir != "":
            self.wDir_file.config(bg="green", activebackground="green", relief=tk.SUNKEN)
            self.output_text.insert(tk.END, f"Working directory set to {self.wDir}\n")
        else:
            self.wDir = None

    def browse_gene(self):
        self.geneFile = filedialog.askopenfilename(title="Select a File")
        if self.geneFile != "":
            self.gene_file.config(bg="green", activebackground="green", relief=tk.SUNKEN)
            self.output_text.insert(tk.END, f"Gene file set to {self.geneFile}\n")
        else:
            self.geneFile = None

    # create a function to input custom mutations
    def custom_mutations_input(self):
        # create a new window
        newWindow = tk.Toplevel(self)
        newWindow.title("Custom Mutations          (close window when finished)")
        newWindow.geometry("700x1000")
        # create a text box for custom mutations
        custom_mutations_window = tk.Text(newWindow, height=50, width=80)
        custom_mutations_window.pack()
        # example mutations
        if self.custom_mutations == {}:
            custom_mutations_window.insert(
                tk.END,
                "Positions:Mutations\n1-10:All\n11-20:A,C,D,E,F,G,H,I,K,L\n33:A,C",
            )
        else:
            custom_mutations_window.insert(tk.END, "Positions:Mutations\n")
            for key, value in self.custom_mutations.items():
                custom_mutations_window.insert(tk.END, str(key) + ":" + value + "\n")

        # create a button to set custom mutations
        def on_closing():
            # set custom mutations
            mutation_text = custom_mutations_window.get("1.0", "end-1c").strip().split("\n")
            self.custom_mutations = parse_custom_mutations(mutation_text[1:])
            newWindow.destroy()

        newWindow.protocol("WM_DELETE_WINDOW", on_closing)
        self.custom_mutations_button.config(bg="green", activebackground="green", relief=tk.SUNKEN)
        self.include_substitutions.set(1)


if __name__ == "__main__":
    root = tk.Tk()
    root.geometry("700x1100")
    app = Application(master=root)
    app.mainloop()
