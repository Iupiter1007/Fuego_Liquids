# mojave_viewer.py
import importlib.util
import os
import tkinter as tk
from tkinter import filedialog, messagebox
from math import pi
import numpy as np
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

# --- Defaults (meters) ---
DEFAULT_GEOMETRY = {
    "ox_d": 0.065,         # oxidizer tank diameter
    "ox_L": 0.30,          # oxidizer tank length
    "fuel_d": 0.065,       # fuel tank diameter
    "fuel_L": 0.30,        # fuel tank length
    "chamber_d": 0.050,    # combustion chamber diameter
    "chamber_L": 0.08,     # combustion chamber length
    "nozzle_throat_d": 0.020,
    "nozzle_exit_d": 0.050,
    "injector_d": 0.020,
    "stack_gap": 0.02      # distance between tanks/chamber
}

class MojaveViewerApp:
    def __init__(self, master):
        self.master = master
        master.title("Mojave Sphinx cross-section viewer")

        # geometry state
        self.geometry = DEFAULT_GEOMETRY.copy()

        # left: controls
        ctrl = tk.Frame(master)
        ctrl.pack(side=tk.LEFT, fill=tk.Y, padx=8, pady=8)

        # import button
        btn_import = tk.Button(ctrl, text="Import .py (your file)", command=self.import_user_py)
        btn_import.pack(fill=tk.X, pady=(0,6))

        # load defaults / reset
        btn_reset = tk.Button(ctrl, text="Reset defaults", command=self.reset_defaults)
        btn_reset.pack(fill=tk.X, pady=(0,6))

        # parameter entries
        self.entries = {}
        for key in ["ox_d","ox_L","fuel_d","fuel_L","chamber_d","chamber_L","nozzle_throat_d","nozzle_exit_d","injector_d","stack_gap"]:
            frm = tk.Frame(ctrl)
            frm.pack(fill=tk.X, pady=2)
            label = tk.Label(frm, text=key.replace("_"," "), width=14, anchor='w')
            label.pack(side=tk.LEFT)
            var = tk.StringVar(value=str(self.geometry[key]))
            ent = tk.Entry(frm, textvariable=var, width=10)
            ent.pack(side=tk.LEFT)
            self.entries[key] = var

        # redraw button
        btn_redraw = tk.Button(ctrl, text="Redraw", command=self.read_inputs_and_draw)
        btn_redraw.pack(fill=tk.X, pady=(6,2))

        # export (png)
        btn_save = tk.Button(ctrl, text="Save PNG", command=self.save_png)
        btn_save.pack(fill=tk.X, pady=(2,6))

        # right: drawing canvas
        fig = Figure(figsize=(6,8))
        self.ax = fig.add_subplot(111)
        self.ax.set_aspect('equal')
        self.ax.axis('off')

        canvas = FigureCanvasTkAgg(fig, master=master)
        canvas.get_tk_widget().pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)
        self.canvas = canvas

        # initial draw
        self.draw_cross_section()

    def reset_defaults(self):
        self.geometry = DEFAULT_GEOMETRY.copy()
        for k,v in self.geometry.items():
            self.entries[k].set(str(v))
        self.draw_cross_section()

    def import_user_py(self):
        path = filedialog.askopenfilename(title="Select your .py file", filetypes=[("Python files","*.py")])
        if not path:
            return
        try:
            geom = self.load_geometry_from_py(path)
            if not geom:
                messagebox.showinfo("Import", f"No GEOMETRY or get_geometry() found in {os.path.basename(path)}.\nYou can still edit parameters manually.")
                return
            # update geometry
            for k in DEFAULT_GEOMETRY.keys():
                if k in geom:
                    self.geometry[k] = float(geom[k])
                    self.entries[k].set(str(self.geometry[k]))
            self.draw_cross_section()
            messagebox.showinfo("Import", f"Loaded geometry from {os.path.basename(path)}")
        except Exception as e:
            messagebox.showerror("Import error", f"Failed to import file:\n{e}")

    def load_geometry_from_py(self, path):
        """
        Attempts to import a module from path and read GEOMETRY or get_geometry()
        """
        spec = importlib.util.spec_from_file_location("user_module_for_geom", path)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        geom = {}
        if hasattr(module, "GEOMETRY"):
            geom = getattr(module, "GEOMETRY")
        elif hasattr(module, "get_geometry"):
            func = getattr(module, "get_geometry")
            geom = func()
        # Validate: keep only numeric keys we expect
        geom_out = {}
        for k in DEFAULT_GEOMETRY.keys():
            if k in geom:
                try:
                    geom_out[k] = float(geom[k])
                except Exception:
                    pass
        return geom_out

    def read_inputs_and_draw(self):
        # read entries into geometry dict
        for k,var in self.entries.items():
            try:
                self.geometry[k] = float(var.get())
            except ValueError:
                messagebox.showerror("Input error", f"Invalid numeric value for {k}: {var.get()}")
                return
        self.draw_cross_section()

    def draw_rect(self, x0, y0, w, h, **kwargs):
        self.ax.add_patch(matplotlib.patches.Rectangle((x0,y0), w, h, **kwargs))

    def draw_cross_section(self):
        self.ax.clear()
        self.ax.set_aspect('equal')
        self.ax.axis('off')

        g = self.geometry
        # stack from top to bottom: oxidizer tank, gap, fuel tank, gap, chamber/nozzle
        # We'll draw centerline at x=0; horizontal cross-section represented as rectangles along y-axis stack.

        # convert diameters to widths for drawing (scale later)
        items = [
            ("Oxidizer tank", g["ox_d"], g["ox_L"]),
            ("gap1", 0.0, g["stack_gap"]),
            ("Fuel tank", g["fuel_d"], g["fuel_L"]),
            ("gap2", 0.0, g["stack_gap"]),
            ("Chamber", g["chamber_d"], g["chamber_L"]),
            ("Nozzle throat", g["nozzle_throat_d"], g["chamber_L"]*0.2),
            ("Nozzle exit", g["nozzle_exit_d"], g["chamber_L"]*0.4),
        ]

        # compute total length
        total_length = sum([h for (_name,d,h) in items])
        # set a drawing scale so that figure fits nicely
        fig_height = total_length
        # x-size: use max diameter * 1.6
        max_diam = max([d for (_name,d,_h) in items] + [0.05])
        # margins
        x_margin = max_diam * 0.6
        y = 0.0
        # We'll draw rectangles centered on x=0 horizontally, with width = diameter
        for name, diam, h in items:
            if diam <= 0:
                # draw empty gap as line
                y += h
                continue
            w = diam
            x0 = -w/2.0
            # draw outer body
            rect = matplotlib.patches.Rectangle((x0,y), w, h, edgecolor='black', facecolor='#ddddff', linewidth=1.2)
            self.ax.add_patch(rect)
            # label
            cx = 0
            cy = y + h/2.0
            self.ax.text(cx, cy, name, ha='center', va='center', fontsize=9)
            y += h

        # draw injector (as small circle) placed at interface between fuel tank and chamber
        # compute y for injector: between end of fuel tank and start of chamber
        y_acc = 0.0
        injector_y = None
        for name, diam, h in items:
            if name == "Fuel tank":
                y_acc += h
            elif name == "gap2":
                injector_y = y_acc + h/2.0
                y_acc += h
            else:
                y_acc += h
        if injector_y is None:
            injector_y = g["ox_L"] + g["stack_gap"] + g["fuel_L"] + 0.5*g["stack_gap"]

        # injector radius roughly injector_d/2
        inj_r = g["injector_d"]/2.0
        # scale injector to be visible even if small
        inj_r_plot = max(inj_r, 0.006)
        circ = matplotlib.patches.Circle((0.0, injector_y), inj_r_plot, edgecolor='red', facecolor='#ffcccc', linewidth=1.0)
        self.ax.add_patch(circ)
        self.ax.text(0.0, injector_y, "Injector", ha='center', va='center', fontsize=8, color='darkred')

        # centerline
        self.ax.plot([-max_diam, max_diam], [0,0], alpha=0.0)  # no visible line, just to set extents

        # autoscale and set limits with margin
        self.ax.set_xlim(-max_diam - x_margin, max_diam + x_margin)
        self.ax.set_ylim(0, total_length + 0.02)
        self.canvas.draw()

    def save_png(self):
        f = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG image","*.png")])
        if not f:
            return
        self.canvas.figure.savefig(f, dpi=200)
        messagebox.showinfo("Saved", f"Saved image to {f}")


if __name__ == "__main__":
    root = tk.Tk()
    app = MojaveViewerApp(root)
    root.mainloop()
