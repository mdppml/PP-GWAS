import os, sys, subprocess, tkinter as tk
from tkinter import ttk, messagebox
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from data_generation import run_all

def resource_path(p):
    if getattr(sys, "frozen", False):
        base = sys._MEIPASS
    else:
        base = os.path.dirname(os.path.abspath(__file__))
    return os.path.join(base, p)

class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("PP-GWAS One-Click")
        self.geometry("900x600")
        self.n_var = tk.StringVar(value="10000")
        self.m_var = tk.StringVar(value="15000")
        self.c_var = tk.StringVar(value="10")
        self.p_var = tk.StringVar(value="4")
        self.b_var = tk.StringVar(value="2")
        self.bpr_var = tk.StringVar(value="")
        self.status_var = tk.StringVar(value="Idle")
        self.canvas = None
        self.build_ui()

    def build_ui(self):
        frm = ttk.Frame(self, padding=12)
        frm.pack(fill="x")
        for i, (lbl, var) in enumerate([("N", self.n_var), ("M", self.m_var), ("C", self.c_var), ("P", self.p_var), ("B", self.b_var)]):
            ttk.Label(frm, text=lbl, width=4).grid(row=0, column=2*i, sticky="w")
            ttk.Entry(frm, textvariable=var, width=10).grid(row=0, column=2*i+1, padx=(0,12))
        ttk.Label(frm, text="BPR", width=4).grid(row=1, column=0, sticky="w", pady=(8,0))
        ttk.Entry(frm, textvariable=self.bpr_var, width=10).grid(row=1, column=1, padx=(0,12), pady=(8,0))
        ttk.Button(frm, text="Run", command=self.run_pipeline).grid(row=1, column=2, pady=(8,0))
        ttk.Label(frm, textvariable=self.status_var).grid(row=1, column=3, columnspan=6, sticky="w", pady=(8,0))
        ttk.Label(frm, text="BPR = blocks per run; defaults to B if left empty. Lower it if you hit memory limits (min 1).", foreground="gray").grid(row=2, column=0, columnspan=9, sticky="w", pady=(4,0))
        self.plot_frame = ttk.Frame(self, padding=12)
        self.plot_frame.pack(fill="both", expand=True)

    def run_pipeline(self):
        try:
            N = int(self.n_var.get()); M = int(self.m_var.get()); C = int(self.c_var.get()); P = int(self.p_var.get()); B = int(self.b_var.get())
        except:
            messagebox.showerror("Error", "Inputs must be integers.")
            return
        BPR = self.bpr_var.get().strip()
        if BPR == "":
            BPR = B
        try:
            BPR = int(BPR)
        except:
            messagebox.showerror("Error", "BPR must be an integer.")
            return
        self.status_var.set("Generating synthetic data...")
        self.update_idletasks()
        try:
            run_all(N, M, C, P, B)
        except Exception as e:
            messagebox.showerror("Error", f"Data generation failed:\n{e}")
            return
        self.status_var.set("Running PP-GWAS...")
        self.update_idletasks()
        script = resource_path("ppgwas.sh")
        if not os.path.exists(script):
            messagebox.showerror("Error", "ppgwas.sh not found.")
            return
        cmd = ["bash", script, "8000", str(N), str(M), str(C), str(P), str(B), "5", str(BPR)]
        try:
            subprocess.run(cmd, check=True, cwd=resource_path("."))
        except Exception as e:
            messagebox.showerror("Error", f"ppgwas.sh failed:\n{e}. Please refer to the logs for more information.")
            return
        self.status_var.set("Loading results...")
        self.update_idletasks()
        file_path = os.path.join(resource_path("."), "Data", f"N{N}_M{M}_C{C}_P{P}_B{B}", "neg_log_transfer.npy")
        if not os.path.exists(file_path):
            messagebox.showerror("Error", f"Result file not found:\n{file_path}")
            return
        try:
            neg_log_p = np.load(file_path)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to load results:\n{e}")
            return
        positions = np.arange(len(neg_log_p))
        fig = plt.Figure(figsize=(10,5), dpi=100)
        ax = fig.add_subplot(111)
        ax.scatter(positions, neg_log_p, s=1)
        ax.set_xlim(0, M)
        ax.set_ylim(bottom=0)
        ax.set_xlabel("SNP Index")
        ax.set_ylabel("-log10(p-value)")
        ax.set_title("Manhattan Plot")
        if self.canvas:
            self.canvas.get_tk_widget().destroy()
        self.canvas = FigureCanvasTkAgg(fig, master=self.plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.status_var.set("Done.")

if __name__ == "__main__":
    App().mainloop()

