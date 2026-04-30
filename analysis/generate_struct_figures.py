#!/usr/bin/env python3
"""Generate publication-ready structural figures for comprehensive_guide.md.

Unified color scheme:
  - Protein: gray90 cartoon
  - L1 (active-site MVAP): marine blue
  - L2 (inhibitory MVAP): firebrick red
  - Entrance residues: tv_yellow
  - Catalytic residues: forest green
  - Q212: orange
  - H-bonds: yellow dashes
"""

import sys, os, tempfile, shutil
sys.argv = ['pymol', '-cq']
import pymol
pymol.finish_launching(['pymol', '-cq'])
from pymol import cmd, stored

ANALYSIS_DIR = "/nfs/lambda_stor_01/homes/rzhu/PMD/analysis"
BOLTZ_2MVAP = "/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/output_2mvap"
BOLTZ_BEST = "/nfs/lambda_stor_01/homes/rzhu/PMD/structures/boltz2/best_structures"

VARIANT_DIRS = {
    'WT': 'PMDsc_WT_2mvap',
    'Y19H': 'PMDsc_Y19H_2mvap',
    'R74G': 'PMDsc_R74G_2mvap',
    'GKQ': 'PMDsc_R74G_R147K_M212Q_2mvap',
    'HKQ': 'PMDsc_R74H_R147K_M212Q_2mvap',
}

WIDTH, HEIGHT = 2400, 1800
DPI = 300


def fix_boltz2_pdb(pdb_path, output_path):
    """Fix Boltz-2 PDB multi-char chain IDs (L2 -> B) and shifted columns."""
    with open(pdb_path) as f:
        lines = f.readlines()
    fixed = []
    for line in lines:
        if not line.startswith(('ATOM', 'HETATM')):
            if line.startswith('TER') and 'L2' in line:
                line = line.replace('L2', 'B ')
            fixed.append(line)
            continue
        raw = line.rstrip('\n')
        try:
            float(raw[30:38]); float(raw[38:46]); float(raw[46:54])
            fixed.append(line)
            continue
        except (ValueError, IndexError):
            pass
        try:
            x = float(raw[31:39]); y = float(raw[39:47]); z = float(raw[47:55])
            record = raw[:6]; serial = raw[6:11]; atom_name = raw[12:16]
            resn = raw[17:20]; resi = raw[24:27].strip() or '1'
            occ = raw[55:61] if len(raw) > 60 else ' 1.00'
            bfac = raw[61:67] if len(raw) > 66 else ' 0.00'
            elem = raw[76:78].strip() if len(raw) > 77 else atom_name.strip()[0]
            new_line = (f"{record}{serial} {atom_name} {resn} B{resi:>4s}    "
                        f"{x:8.3f}{y:8.3f}{z:8.3f}{occ}{bfac}          {elem:>2s}  \n")
            fixed.append(new_line)
        except (ValueError, IndexError):
            fixed.append(line)
    with open(output_path, 'w') as f:
        f.writelines(fixed)


def get_2mvap_pdb(variant, model_idx):
    vdir = VARIANT_DIRS[variant]
    return (f"{BOLTZ_2MVAP}/{vdir}/boltz_results_{vdir}/"
            f"predictions/{vdir}/{vdir}_model_{model_idx}.pdb")


def get_best_pdb(variant):
    name_map = {
        'WT': 'WT', 'Y19H': 'Y19H', 'R74G': 'R74G',
        'GKQ': 'R74G_R147K_M212Q', 'HKQ': 'R74H_R147K_M212Q',
    }
    return f"{BOLTZ_BEST}/{name_map[variant]}.pdb"


def setup_common():
    cmd.set("ray_opaque_background", 1)
    cmd.set("bg_rgb", [1.0, 1.0, 1.0])
    cmd.set("ray_trace_mode", 0)
    cmd.set("antialias", 2)
    cmd.set("ray_shadows", 0)
    cmd.set("depth_cue", 0)
    cmd.set("fog", 0)
    cmd.set("specular", 0.3)
    cmd.set("spec_reflect", 0.3)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("cartoon_flat_sheets", 1)
    cmd.set("stick_radius", 0.2)
    cmd.set("sphere_scale", 0.3)
    cmd.set("label_size", 11)
    cmd.set("label_font_id", 7)
    cmd.set("label_color", "black")
    cmd.set("label_outline_color", "white")
    cmd.set("label_position", [0, 1, 3])
    cmd.set("dash_width", 3.0)
    cmd.set("dash_gap", 0.2)
    cmd.set("dash_length", 0.15)
    cmd.set("dash_color", "tv_yellow")


def label_res(sel, text):
    cmd.label(f"({sel}) and name CA", f'"{text}"')


def label_lig(sel, text):
    stored.idx = None
    cmd.iterate(f"({sel}) and name P1", "stored.idx = index", quiet=1)
    if stored.idx is None:
        cmd.iterate(f"({sel}) and elem P", "stored.idx = index", quiet=1)
    if stored.idx is None:
        cmd.iterate(f"({sel})", "stored.idx = index", quiet=1)
    if stored.idx:
        cmd.label(f"index {stored.idx}", f'"{text}"')


def select_l2(obj_name):
    cmd.select("L2", f"{obj_name} and chain B and resn LIG")
    if cmd.count_atoms("L2") == 0:
        cmd.select("L2", f"{obj_name} and chain L and resn LIG and resi 21")


def show_residue_sticks(sel, color, radius=0.2):
    """Show sidechain + CA sticks for a residue selection, colored by element."""
    cmd.show("sticks", f"({sel}) and (sidechain or name CA)")
    cmd.color(color, f"({sel}) and elem C")
    cmd.color("tv_red", f"({sel}) and elem O")
    cmd.color("tv_blue", f"({sel}) and elem N")
    cmd.set("stick_radius", radius, sel)


def show_ligand(sel, carbon_color, radius=0.18, sphere_scale=0.2):
    """Show ligand as ball-and-stick."""
    cmd.show("sticks", sel)
    cmd.show("spheres", sel)
    cmd.color(carbon_color, f"({sel}) and elem C")
    cmd.color("tv_red", f"({sel}) and elem O")
    cmd.color("tv_blue", f"({sel}) and elem N")
    cmd.color("orange", f"({sel}) and elem P")
    cmd.set("stick_radius", radius, sel)
    cmd.set("sphere_scale", sphere_scale, sel)


# ============================================================================
# Figure 1: Overview - WT with L2 at entrance
# ============================================================================
def fig_overview_entrance():
    print("Generating struct_overview_entrance.png...")
    cmd.delete("all")
    setup_common()

    tmpdir = tempfile.mkdtemp()
    pdb_fixed = os.path.join(tmpdir, "wt_fixed.pdb")
    fix_boltz2_pdb(get_2mvap_pdb('WT', 0), pdb_fixed)
    cmd.load(pdb_fixed, "wt")

    cmd.select("prot", "wt and chain A")
    cmd.select("L1", "wt and chain L and resn LIG and resi 1")
    select_l2("wt")
    cmd.select("atp", "wt and chain T")
    cmd.select("mg", "wt and chain M")

    entrance_resi = "22+28+72+74+209"
    cmd.select("entrance_res", f"prot and resi {entrance_resi}")

    cmd.hide("everything")

    # Solid protein cartoon with slight transparency
    cmd.show("cartoon", "prot")
    cmd.color("gray90", "prot")
    cmd.set("cartoon_transparency", 0.15, "prot")

    # Color entrance residues yellow in cartoon
    cmd.color("tv_yellow", f"prot and resi {entrance_resi}")

    # R74 as sticks for emphasis
    show_residue_sticks("prot and resi 74", "tv_yellow", 0.25)

    # L1 inside the pocket: blue ball-and-stick
    show_ligand("L1", "marine", radius=0.22, sphere_scale=0.28)
    # L2 at the entrance: red ball-and-stick
    show_ligand("L2", "firebrick", radius=0.22, sphere_scale=0.28)

    # ATP thin sticks for context
    cmd.show("sticks", "atp")
    cmd.color("palecyan", "atp")
    cmd.set("stick_radius", 0.08, "atp")

    # Mg sphere
    cmd.show("spheres", "mg")
    cmd.color("green", "mg")
    cmd.set("sphere_scale", 0.5, "mg")

    # Orient on entrance, looking into the pocket
    cmd.orient("entrance_res")
    cmd.turn("x", -10)
    cmd.zoom("entrance_res or L1 or L2", buffer=6)
    cmd.clip("slab", 40)

    # Labels well separated — L1 above-left, L2 below-right
    cmd.set("label_position", [-5, 5, 5])
    label_lig("L1", "L1")
    cmd.set("label_position", [5, -5, 5])
    label_lig("L2", "L2")
    cmd.set("label_position", [4, 3, 4])
    label_res("prot and resi 74", "R74")

    cmd.deselect()
    cmd.ray(WIDTH, HEIGHT)
    cmd.png(f"{ANALYSIS_DIR}/struct_overview_entrance.png", dpi=DPI)
    shutil.rmtree(tmpdir)
    print("  Done.")


# ============================================================================
# Figure 2: R74 Gateway - WT vs HKQ superposed, focus on loop 17-33
# ============================================================================
def fig_r74_gateway():
    print("Generating struct_r74_gateway.png...")
    cmd.delete("all")
    setup_common()

    cmd.load(get_best_pdb('WT'), "wt")
    cmd.load(get_best_pdb('HKQ'), "hkq")
    cmd.align("hkq", "wt")

    cmd.hide("everything")

    cmd.show("cartoon", "wt and chain A")
    cmd.color("palecyan", "wt and chain A")
    cmd.set("cartoon_transparency", 0.35, "wt and chain A")
    cmd.color("deepteal", "wt and chain A and resi 17-33")

    cmd.show("cartoon", "hkq and chain A")
    cmd.color("lightorange", "hkq and chain A")
    cmd.set("cartoon_transparency", 0.35, "hkq and chain A")
    cmd.color("salmon", "hkq and chain A and resi 17-33")

    show_residue_sticks("wt and chain A and resi 74", "deepteal", 0.28)
    show_residue_sticks("hkq and chain A and resi 74", "salmon", 0.28)

    # Label at sidechain tips to separate spatially
    # R74 (Arg): label at NH2 — tip of the long arginine sidechain
    cmd.set("label_position", [2, 3, 5])
    stored.idx = None
    cmd.iterate("wt and chain A and resi 74 and name NH2", "stored.idx = index", quiet=1)
    if not stored.idx:
        cmd.iterate("wt and chain A and resi 74 and name NH1", "stored.idx = index", quiet=1)
    if stored.idx:
        cmd.label(f"index {stored.idx}", '"R74 (WT)"')
    # H74 (His): label at CA — short sidechain, label near backbone
    cmd.set("label_position", [-3, -4, 5])
    label_res("hkq and chain A and resi 74", "H74 (HKQ)")

    cmd.orient("(wt or hkq) and chain A and resi 17-33+74")
    cmd.turn("y", 25)
    cmd.turn("x", -10)
    cmd.zoom("(wt or hkq) and chain A and resi 17-33+74", buffer=4)
    cmd.clip("slab", 45)

    cmd.deselect()
    cmd.ray(WIDTH, HEIGHT)
    cmd.png(f"{ANALYSIS_DIR}/struct_r74_gateway.png", dpi=DPI)
    print("  Done.")


# ============================================================================
# Figure 3: R74G Competitive binding - L2 inside the active site
# ============================================================================
def fig_r74g_competitive():
    print("Generating struct_r74g_competitive.png...")
    cmd.delete("all")
    setup_common()

    tmpdir = tempfile.mkdtemp()
    pdb_fixed = os.path.join(tmpdir, "r74g_fixed.pdb")
    fix_boltz2_pdb(get_2mvap_pdb('R74G', 1), pdb_fixed)
    cmd.load(pdb_fixed, "r74g")

    cmd.select("prot", "r74g and chain A")
    cmd.select("L1", "r74g and chain L and resn LIG and resi 1")
    select_l2("r74g")
    cmd.select("atp", "r74g and chain T")

    cmd.hide("everything")

    # Protein cartoon — visible context
    cmd.show("cartoon", "prot")
    cmd.color("gray90", "prot")
    cmd.set("cartoon_transparency", 0.15, "prot")

    # G74: yellow sticks (no sidechain — just backbone)
    show_residue_sticks("prot and resi 74", "tv_yellow", 0.2)
    # Y19: green sticks (catalytic context)
    show_residue_sticks("prot and resi 19", "forest", 0.2)

    # L1: marine blue
    show_ligand("L1", "marine", radius=0.22, sphere_scale=0.28)
    # L2: firebrick red
    show_ligand("L2", "firebrick", radius=0.22, sphere_scale=0.28)

    # ATP thin sticks for context
    cmd.show("sticks", "atp")
    cmd.color("palecyan", "atp")
    cmd.set("stick_radius", 0.08, "atp")

    # Labels well separated
    cmd.set("label_position", [-5, 4, 5])
    label_lig("L1", "L1")
    cmd.set("label_position", [5, -4, 5])
    label_lig("L2", "L2")
    cmd.set("label_position", [4, 3, 4])
    label_res("prot and resi 74", "G74")
    cmd.set("label_position", [-4, -2, 4])
    label_res("prot and resi 19", "Y19")

    cmd.orient("L1 or L2")
    cmd.turn("y", 15)
    cmd.zoom("L1 or L2", buffer=6)
    cmd.clip("slab", 30)

    cmd.deselect()
    cmd.ray(WIDTH, HEIGHT)
    cmd.png(f"{ANALYSIS_DIR}/struct_r74g_competitive.png", dpi=DPI)
    shutil.rmtree(tmpdir)
    print("  Done.")


# ============================================================================
# Figure 4: Q212 bridge in HKQ
# ============================================================================
def fig_q212_bridge():
    print("Generating struct_q212_bridge.png...")
    cmd.delete("all")
    setup_common()

    cmd.load(get_best_pdb('HKQ'), "hkq")
    cmd.select("prot", "hkq and chain A")

    cmd.hide("everything")
    cmd.show("cartoon", "prot")
    cmd.color("gray90", "prot")
    cmd.set("cartoon_transparency", 0.15, "prot")

    # Q212 bridge residues
    show_residue_sticks("prot and resi 212", "orange", 0.25)
    show_residue_sticks("prot and resi 19", "forest", 0.22)
    show_residue_sticks("prot and resi 22", "forest", 0.22)
    show_residue_sticks("prot and resi 74", "tv_yellow", 0.2)

    # K147-E144 salt bridge residues
    show_residue_sticks("prot and resi 147", "purple", 0.22)
    show_residue_sticks("prot and resi 144", "purple", 0.22)

    # Q212 hydrogen bonds (yellow dashes)
    cmd.distance("q212_y19",
                 "prot and resi 212 and (name OE1 or name NE2)",
                 "prot and resi 19 and (name OH or name CE2)",
                 mode=0)
    cmd.color("tv_yellow", "q212_y19")
    cmd.hide("labels", "q212_y19")

    cmd.distance("q212_k22",
                 "prot and resi 212 and (name OE1 or name NE2)",
                 "prot and resi 22 and name NZ",
                 mode=0)
    cmd.color("tv_yellow", "q212_k22")
    cmd.hide("labels", "q212_k22")

    # K147-E144 salt bridge (magenta dashes)
    cmd.distance("k147_e144",
                 "prot and resi 147 and name NZ",
                 "prot and resi 144 and (name OE1 or name OE2)",
                 mode=0)
    cmd.color("magenta", "k147_e144")
    cmd.hide("labels", "k147_e144")

    # Labels
    cmd.set("label_position", [0, 4, 4])
    label_res("prot and resi 212", "Q212")
    cmd.set("label_position", [-3, -2, 3])
    label_res("prot and resi 19", "Y19")
    cmd.set("label_position", [3, -2, 3])
    label_res("prot and resi 22", "K22")
    cmd.set("label_position", [4, 3, 3])
    label_res("prot and resi 74", "H74")
    cmd.set("label_position", [-3, 3, 3])
    label_res("prot and resi 147", "K147")
    cmd.set("label_position", [-4, -1, 3])
    label_res("prot and resi 144", "E144")

    # Widen view to include K147/E144
    cmd.orient("prot and resi 19+22+74+144+147+212")
    cmd.turn("y", -20)
    cmd.turn("x", 10)
    cmd.zoom("prot and resi 19+22+74+144+147+212", buffer=3)
    cmd.clip("slab", 30)

    cmd.deselect()
    cmd.ray(WIDTH, HEIGHT)
    cmd.png(f"{ANALYSIS_DIR}/struct_q212_bridge.png", dpi=DPI)
    print("  Done.")


if __name__ == "__main__":
    fig_overview_entrance()
    fig_r74_gateway()
    fig_r74g_competitive()
    fig_q212_bridge()
    print(f"\nAll 4 figures generated in {ANALYSIS_DIR}/")
    cmd.quit()
