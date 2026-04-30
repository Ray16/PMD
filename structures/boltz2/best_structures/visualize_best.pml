# PyMOL script: PMDsc Boltz-2 cofolded structures — catalytic residues highlighted
# Usage: cd to best_structures/ then run: pymol visualize_best.pml

# --- Load best structures ---
load WT.pdb, WT
load Y19H.pdb, Y19H
load R74G.pdb, R74G
load R74H_R147K_M212Q.pdb, HKQ
load R74G_R147K_M212Q.pdb, GKQ

# ============================================================================
# Rendering settings (dark background)
# ============================================================================
bg_color black
set antialias, 2
set ambient, 0.45
set direct, 0.5
set reflect, 0.4
set specular, 0.4
set shininess, 25
set ray_shadow, 0
set ray_trace_mode, 0
set depth_cue, 1
set fog_start, 0.55
set fog, 0.3

# Cartoon settings
set cartoon_smooth_loops, 1
set cartoon_flat_sheets, 1
set cartoon_fancy_helices, 1
set cartoon_highlight_color, grey70
set cartoon_side_chain_helper, 1
set cartoon_oval_length, 1.2
set cartoon_loop_radius, 0.3
set cartoon_transparency, 0.15

# Stick settings
set stick_radius, 0.15
set stick_ball, 0
set valence, 0

# Label settings
set label_font_id, 7
set label_size, 14
set label_outline_color, black
set label_color, white
set label_position, [0, 0, 3]

# ============================================================================
# Protein display
# ============================================================================
hide everything, all
show cartoon, chain A
color grey50, chain A

# ============================================================================
# Catalytic residues (11 verified, colored by functional role)
# ============================================================================

# Catalytic core: D302 (catalytic base), R158 (decarboxylation)
select cat_core, resi 302+158 and chain A
show sticks, cat_core
util.cbam cat_core

# ATP binding: S121, S155, S208
select atp_site, resi 121+155+208 and chain A
show sticks, atp_site
util.cbao atp_site

# Substrate binding: K18, Y19, W20, K22, S120, S153
select sub_site, resi 18+19+20+22+120+153 and chain A
show sticks, sub_site
util.cbac sub_site

# ============================================================================
# Ligands
# ============================================================================

# MVAP (chain L) — green carbons
show sticks, chain L
util.cbag chain L
set stick_radius, 0.18, chain L

# ATP (chain T) — pink carbons
show sticks, chain T
util.cbap chain T
set stick_radius, 0.18, chain T

# Mg2+ (chain M) — bright sphere
show spheres, chain M
color limegreen, chain M
set sphere_scale, 0.6, chain M

# ============================================================================
# Labels on WT only
# ============================================================================
label WT and chain A and resi 302 and name CA, "D302"
label WT and chain A and resi 158 and name CA, "R158"
label WT and chain A and resi 121 and name CA, "S121"
label WT and chain A and resi 155 and name CA, "S155"
label WT and chain A and resi 208 and name CA, "S208"
label WT and chain A and resi 18 and name CA, "K18"
label WT and chain A and resi 19 and name CA, "Y19"
label WT and chain A and resi 20 and name CA, "W20"
label WT and chain A and resi 22 and name CA, "K22"
label WT and chain A and resi 120 and name CA, "S120"
label WT and chain A and resi 153 and name CA, "S153"

# ============================================================================
# View setup
# ============================================================================
disable Y19H
disable R74G
disable HKQ
disable GKQ

zoom resi 18+19+20+22+120+121+153+155+158+208+302 and chain A and WT, 8

# Deselect all to clean up
deselect

# ============================================================================
# Render high-resolution image:
#   ray 2400, 2400
#   png PMDsc_WT_active_site.png, dpi=300
#
# Superimpose variants:
#   align Y19H, WT
#   align R74G, WT
#   align HKQ, WT
#   align GKQ, WT
#
# Show mutation sites (not shown by default):
#   show sticks, resi 74+147+212 and chain A
#   util.cbay resi 74+147+212 and chain A
#
# Color key (CPK heteroatoms: red=O, blue=N, orange=P, yellow=S):
#   Magenta carbons:  Catalytic core (D302, R158)
#   Orange carbons:   ATP binding (S121, S155, S208)
#   Cyan carbons:     Substrate binding (K18, Y19, W20, K22, S120, S153)
#   Green carbons:    MVAP ligand (chain L)
#   Pink carbons:     ATP ligand (chain T)
#   Green sphere:     Mg2+ ion
#   Grey:             Protein cartoon
# ============================================================================
