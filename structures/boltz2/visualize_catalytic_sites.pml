# PyMOL script: Visualize PMDsc Boltz-2 cofolding results with catalytic residues
# Usage: pymol visualize_catalytic_sites.pml

# --- Load all 5 structures ---
load output_1sample/PMDsc_WT/boltz_results_PMDsc_WT/predictions/PMDsc_WT/PMDsc_WT_model_0.pdb, WT
load output_1sample/PMDsc_Y19H/boltz_results_PMDsc_Y19H/predictions/PMDsc_Y19H/PMDsc_Y19H_model_0.pdb, Y19H
load output_1sample/PMDsc_R74G/boltz_results_PMDsc_R74G/predictions/PMDsc_R74G/PMDsc_R74G_model_0.pdb, R74G
load output_1sample/PMDsc_R74H_R147K_M212Q/boltz_results_PMDsc_R74H_R147K_M212Q/predictions/PMDsc_R74H_R147K_M212Q/PMDsc_R74H_R147K_M212Q_model_0.pdb, HKQ
load output_1sample/PMDsc_R74G_R147K_M212Q/boltz_results_PMDsc_R74G_R147K_M212Q/predictions/PMDsc_R74G_R147K_M212Q/PMDsc_R74G_R147K_M212Q_model_0.pdb, GKQ

# --- General display ---
bg_color white
set ray_shadow, 0
set cartoon_fancy_helices, 1
set cartoon_side_chain_helper, 1

# Show protein as cartoon
hide everything, all
show cartoon, all
color palecyan, WT and chain A
color palegreen, Y19H and chain A
color lightorange, R74G and chain A
color lightpink, HKQ and chain A
color lightblue, GKQ and chain A

# --- Catalytic residues (sticks, colored by role) ---
# Catalytic core: D302, R158
select catalytic_core, (resi 302+158) and chain A
show sticks, catalytic_core
color red, catalytic_core

# ATP binding: S121, S155, S208
select atp_binding, (resi 121+155+208) and chain A
show sticks, atp_binding
color orange, atp_binding

# Substrate binding: K18, Y19, W20, K22, S120, S153
select substrate_binding, (resi 18+19+20+22+120+153) and chain A
show sticks, substrate_binding
color marine, substrate_binding

# --- Mutation sites (thicker sticks, distinct color) ---
select mutation_sites, (resi 19+74+147+212) and chain A
show sticks, mutation_sites
color yellow, mutation_sites
set stick_radius, 0.2, mutation_sites

# --- Ligands ---
# MVAP (chain B typically), ATP (chain C), Mg2+ (chain D)
select mvap, not chain A and not resn MG and not resn ATP
select atp_lig, resn ATP
select mg_ion, resn MG

show sticks, mvap
show sticks, atp_lig
show spheres, mg_ion

color green, mvap
color cyan, atp_lig
color magenta, mg_ion
set sphere_scale, 0.5, mg_ion

# --- Labels for catalytic residues (WT only) ---
set label_size, 12
set label_color, black
label WT and chain A and resi 302 and name CA, "D302 (base)"
label WT and chain A and resi 158 and name CA, "R158 (decarb)"
label WT and chain A and resi 18 and name CA, "K18"
label WT and chain A and resi 19 and name CA, "Y19"
label WT and chain A and resi 20 and name CA, "W20"
label WT and chain A and resi 22 and name CA, "K22"
label WT and chain A and resi 74 and name CA, "R74"
label WT and chain A and resi 121 and name CA, "S121"
label WT and chain A and resi 147 and name CA, "R147"
label WT and chain A and resi 153 and name CA, "S153"
label WT and chain A and resi 155 and name CA, "S155"
label WT and chain A and resi 208 and name CA, "S208"
label WT and chain A and resi 212 and name CA, "M212"

# --- Start with WT visible, others hidden ---
disable Y19H
disable R74G
disable HKQ
disable GKQ

# --- Zoom to active site ---
select active_site, (resi 18+19+20+22+74+120+121+147+153+155+158+208+212+302) and chain A and WT
zoom active_site, 8

# --- Save session ---
# save PMDsc_boltz2_catalytic.pse

# --- Legend ---
# Red sticks:     Catalytic core (D302, R158)
# Orange sticks:  ATP binding (S121, S155, S208)
# Blue sticks:    Substrate binding (K18, Y19, W20, K22, S120, S153)
# Yellow sticks:  Mutation sites (Y19, R74, R147, M212)
# Green sticks:   MVAP ligand
# Cyan sticks:    ATP ligand
# Magenta sphere: Mg2+ ion
#
# Toggle variants on/off in the object panel:
#   WT, Y19H, R74G, HKQ (R74H-R147K-M212Q), GKQ (R74G-R147K-M212Q)
