#!/usr/bin/env python3
"""Generate enzyme overview (no labels, with legend) and 5-variant panel figure."""

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

PANEL_MODELS = {
    'WT': 0,
    'Y19H': 2,
    'R74G': 1,
    'GKQ': 9,
    'HKQ': 0,
}


def fix_boltz2_pdb(pdb_path, output_path):
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
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_smooth_loops", 1)
    cmd.set("cartoon_flat_sheets", 1)
    cmd.set("label_size", 11)
    cmd.set("label_font_id", 7)
    cmd.set("label_color", "black")
    cmd.set("label_outline_color", "white")
    cmd.set("label_position", [0, 1, 3])


# ============================================================================
# Enzyme overview: WT structure, NO text labels, color-coded with legend
# ============================================================================
def fig_enzyme_overview():
    print("Generating struct_enzyme_overview.png...")
    cmd.delete("all")
    setup_common()

    cmd.load(get_best_pdb('WT'), "wt")
    cmd.select("prot", "wt and chain A")

    cmd.hide("everything")
    cmd.show("cartoon", "prot")
    cmd.color("gray90", "prot")

    # Color loop 17-33 in cartoon
    cmd.color("lightorange", "prot and resi 17-33")

    # Catalytic residues: green sticks (sidechain + CA)
    cat_resi = "18+19+20+22+120+121+153+155+158+208+302"
    cmd.show("sticks", f"prot and resi {cat_resi} and (sidechain or name CA)")
    cmd.color("forest", f"prot and resi {cat_resi} and elem C")
    cmd.color("tv_red", f"prot and resi {cat_resi} and elem O")
    cmd.color("tv_blue", f"prot and resi {cat_resi} and elem N")
    cmd.set("stick_radius", 0.2, f"prot and resi {cat_resi}")

    # R74: yellow sticks
    cmd.show("sticks", "prot and resi 74 and (sidechain or name CA)")
    cmd.color("tv_yellow", "prot and resi 74 and elem C")
    cmd.color("tv_red", "prot and resi 74 and elem O")
    cmd.color("tv_blue", "prot and resi 74 and elem N")
    cmd.set("stick_radius", 0.25, "prot and resi 74")

    # R147: purple sticks
    cmd.show("sticks", "prot and resi 147 and (sidechain or name CA)")
    cmd.color("purple", "prot and resi 147 and elem C")
    cmd.color("tv_red", "prot and resi 147 and elem O")
    cmd.color("tv_blue", "prot and resi 147 and elem N")
    cmd.set("stick_radius", 0.25, "prot and resi 147")

    # M212: orange sticks
    cmd.show("sticks", "prot and resi 212 and (sidechain or name CA)")
    cmd.color("orange", "prot and resi 212 and elem C")
    cmd.color("tv_red", "prot and resi 212 and elem O")
    cmd.color("tv_blue", "prot and resi 212 and elem N")
    cmd.set("stick_radius", 0.25, "prot and resi 212")

    # Entrance residues N72, N28, T209: yellow sticks
    for resi in [28, 72, 209]:
        cmd.show("sticks", f"prot and resi {resi} and (sidechain or name CA)")
        cmd.color("tv_yellow", f"prot and resi {resi} and elem C")
        cmd.color("tv_red", f"prot and resi {resi} and elem O")
        cmd.color("tv_blue", f"prot and resi {resi} and elem N")
        cmd.set("stick_radius", 0.2, f"prot and resi {resi}")

    # Orient so entrance (R74, loop 17-33) faces the viewer
    cmd.select("entrance", "prot and resi 22+28+72+74+209")
    cmd.orient("entrance")
    cmd.turn("y", 180)
    cmd.zoom("prot", buffer=3)

    cmd.deselect()

    # Render raw PyMOL image
    tmpdir = tempfile.mkdtemp()
    raw_path = os.path.join(tmpdir, "raw_overview.png")
    cmd.ray(2400, 1800)
    cmd.png(raw_path, dpi=300)

    # Add matplotlib legend
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg
    from matplotlib.patches import Patch

    img = mpimg.imread(raw_path)
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.imshow(img)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axis('off')

    legend_elements = [
        Patch(facecolor='#339933', label='Catalytic residues'),
        Patch(facecolor='#CCCC00', label='Entrance / R74'),
        Patch(facecolor='#FF8000', label='M212 (mutation target)'),
        Patch(facecolor='#BF00BF', label='R147 (mutation target)'),
        Patch(facecolor='#FFCC80', label='Loop 17-33'),
    ]

    fig.legend(handles=legend_elements, loc='lower center', ncol=5,
               fontsize=11, frameon=True, bbox_to_anchor=(0.5, 0.01),
               edgecolor='gray', fancybox=True)

    out_path = f"{ANALYSIS_DIR}/struct_enzyme_overview.png"
    fig.savefig(out_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    shutil.rmtree(tmpdir)
    print("  Done.")


# ============================================================================
# 5-variant panel: one panel per variant showing L2 position (zoomed to entrance)
# ============================================================================
def fig_variant_panels():
    print("Generating 5 variant panel images...")
    tmpdir = tempfile.mkdtemp()

    variants_order = ['WT', 'Y19H', 'R74G', 'GKQ', 'HKQ']
    info = {
        'WT':   {'ki': '18',  'close': '16/50', 'mode': 'entrance'},
        'Y19H': {'ki': 'N.D.', 'close': '8/50', 'mode': 'entrance'},
        'R74G': {'ki': '110', 'close': '9/50', 'mode': 'active site'},
        'GKQ':  {'ki': '11',  'close': '3/50', 'mode': 'entrance'},
        'HKQ':  {'ki': '80',  'close': '1/50', 'mode': 'blocked'},
    }

    cmd.delete("all")
    setup_common()

    wt_fixed = os.path.join(tmpdir, "wt_ref.pdb")
    fix_boltz2_pdb(get_2mvap_pdb('WT', 0), wt_fixed)
    cmd.load(wt_fixed, "wt_ref")

    # Tighter zoom on entrance region
    cmd.select("entrance_ref", "wt_ref and chain A and resi 18+19+22+72+74+209")
    cmd.orient("entrance_ref")
    cmd.turn("y", 20)
    cmd.turn("x", -15)
    cmd.zoom("entrance_ref", buffer=5)
    ref_view = cmd.get_view()

    panel_paths = []
    for variant in variants_order:
        cmd.delete("all")
        setup_common()

        pdb_src = get_2mvap_pdb(variant, PANEL_MODELS[variant])
        pdb_fixed = os.path.join(tmpdir, f"{variant}_fixed.pdb")
        fix_boltz2_pdb(pdb_src, pdb_fixed)
        cmd.load(pdb_fixed, "mol")

        cmd.load(wt_fixed, "ref_align")
        cmd.align("mol and chain A", "ref_align and chain A")
        cmd.delete("ref_align")

        cmd.select("prot", "mol and chain A")
        cmd.select("L1", "mol and chain L and resn LIG and resi 1")
        cmd.select("L2", "mol and chain B and resn LIG")
        if cmd.count_atoms("L2") == 0:
            cmd.select("L2", "mol and chain L and resn LIG and resi 21")
        cmd.select("atp", "mol and chain T")

        cmd.hide("everything")

        cmd.show("cartoon", "prot")
        cmd.color("gray90", "prot")
        cmd.set("cartoon_transparency", 0.2, "prot")

        # Color entrance residues in cartoon only (no sticks)
        entrance_resi = "22+28+72+74+209"
        cmd.color("tv_yellow", f"prot and resi {entrance_resi}")

        cmd.show("sticks", "L1")
        cmd.show("spheres", "L1")
        cmd.color("marine", "L1 and elem C")
        cmd.color("tv_red", "L1 and elem O")
        cmd.color("orange", "L1 and elem P")
        cmd.set("stick_radius", 0.18, "L1")
        cmd.set("sphere_scale", 0.22, "L1")

        cmd.show("sticks", "L2")
        cmd.show("spheres", "L2")
        cmd.color("firebrick", "L2 and elem C")
        cmd.color("tv_red", "L2 and elem O")
        cmd.color("orange", "L2 and elem P")
        cmd.set("stick_radius", 0.18, "L2")
        cmd.set("sphere_scale", 0.22, "L2")

        cmd.show("sticks", "atp")
        cmd.color("palecyan", "atp")
        cmd.set("stick_radius", 0.08, "atp")

        cmd.set_view(ref_view)
        cmd.clip("slab", 30)

        cmd.deselect()
        panel_path = os.path.join(tmpdir, f"panel_{variant}.png")
        cmd.ray(960, 720)
        cmd.png(panel_path, dpi=150)
        panel_paths.append((variant, panel_path))
        print(f"  {variant} panel rendered.")

    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    fig, axes = plt.subplots(1, 5, figsize=(25, 6))
    fig.subplots_adjust(wspace=0.05, left=0.02, right=0.98, top=0.88, bottom=0.02)

    for ax, (variant, path) in zip(axes, panel_paths):
        img = mpimg.imread(path)
        ax.imshow(img)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(1.5)
            spine.set_color('gray')

        v = info[variant]
        ki_str = f"Ki = {v['ki']} mM" if v['ki'] != 'N.D.' else "Ki = N.D."
        title = f"{variant}\n{ki_str} | {v['close']} close | {v['mode']}"
        ax.set_title(title, fontsize=13, fontweight='bold', pad=8)

    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='#0080FF', label='L1 (substrate)'),
        Patch(facecolor='#B22222', label='L2 (inhibitor)'),
        Patch(facecolor='#CCCC00', label='Entrance residues'),
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=3,
               fontsize=12, frameon=True, bbox_to_anchor=(0.5, -0.02))

    out_path = f"{ANALYSIS_DIR}/struct_variant_panels.png"
    fig.savefig(out_path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Combined panel saved to {out_path}")

    shutil.rmtree(tmpdir)
    print("  Done.")


if __name__ == "__main__":
    fig_enzyme_overview()
    fig_variant_panels()
    cmd.quit()
