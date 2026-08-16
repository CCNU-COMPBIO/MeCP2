# ================================================================
# VMD Tcl script
# 3C2I MeCP2 MBD structural classification by SASA
#
# Final classification principle:
#
#   Exposure_ratio = SASA_protein / SASA_isolated
#
#   DNA_buried_ratio =
#       (SASA_protein - SASA_complex) / SASA_protein
#
#   if DNA_buried_ratio > 0.01:
#       Protein-DNA interface
#   elseif Exposure_ratio < 0.20:
#       Protein interior
#   else:
#       Protein surface
#
# Priority:
#   Protein-DNA interface > Protein interior > Protein surface
#
# 3C2I defaults:
#   Protein chain = A
#   DNA chains    = B,C
#   MBD range     = residues 78-162
#
# MSE handling:
#   MSE is treated as MET IN MEMORY before SASA calculation.
#   - resname MSE -> MET
#   - atom name SE -> SD
#   - SE atom radius explicitly set to 1.80 A (sulfur-like)
#   Coordinates are not changed.
#
# No Biopython is required.
#
# Run:
#   vmd -dispdev text -e classify_3C2I_MBD_SASA_VMD.tcl
#
# Or inside VMD Tk Console:
#   source classify_3C2I_MBD_SASA_VMD.tcl
# ================================================================


# -----------------------------
# User settings
# -----------------------------
set pdbfile "3C2I.pdb"

set protein_chain "A"
set dna_chains "B C"

set residue_start 78
set residue_end   162

set probe_radius 1.4
set sasa_samples 960

# Final cutoffs used in the current study
set dna_buried_ratio_cutoff 0.01
set exposure_cutoff          0.20

set output_csv     "3C2I_MBD_VMD_SASA_classification_final.csv"
set output_summary "3C2I_MBD_VMD_SASA_classification_final_summary.txt"


# ================================================================
# Helper procedures
# ================================================================

proc safe_delete {sel} {
    if {$sel ne ""} {
        catch {$sel delete}
    }
}

proc classify_category {dna_buried_ratio exposure dna_cutoff exposure_cutoff} {

    # Priority 1: Protein-DNA interface
    if {$dna_buried_ratio > $dna_cutoff} {
        return "Protein-DNA interface"
    }

    # Priority 2: Protein interior
    if {$exposure < $exposure_cutoff} {
        return "Protein interior"
    }

    # Otherwise: Protein surface
    return "Protein surface"
}


# ================================================================
# Load structure
# ================================================================

if {![file exists $pdbfile]} {
    puts stderr "ERROR: Cannot find $pdbfile"
    puts stderr "Put 3C2I.pdb in the current directory or edit 'pdbfile' in the script."
    exit 1
}

mol new $pdbfile type pdb waitfor all
set molid [molinfo top]

puts "============================================================"
puts "3C2I MeCP2 SASA classification using VMD"
puts "============================================================"
puts "PDB file            : $pdbfile"
puts "Protein chain       : $protein_chain"
puts "DNA chains          : $dna_chains"
puts "Requested range     : $residue_start-$residue_end"
puts "Probe radius        : $probe_radius A"
puts "SASA samples        : $sasa_samples"
puts "Interface metric    : (SASA_protein - SASA_complex) / SASA_protein"
puts "Interface cutoff    : DNA buried ratio > $dna_buried_ratio_cutoff"
puts "Interior cutoff     : Exposure ratio < $exposure_cutoff"
puts ""


# ================================================================
# Treat MSE as MET for SASA calculation
# ================================================================
#
# VMD measure sasa uses the atom radii stored for the molecule.
# Therefore we explicitly change the selenium atom to a sulfur-like
# atom for the SASA calculation.
#
# The original PDB file on disk is NOT modified.
# ================================================================

set mse_all [atomselect $molid "chain $protein_chain and resname MSE"]
set mse_se  [atomselect $molid "chain $protein_chain and resname MSE and name SE"]

set nmse [$mse_all num]

if {$nmse > 0} {

    set mse_resids [lsort -unique [$mse_all get resid]]

    puts "MSE residues found: $mse_resids"
    puts "Treating these residues as MET for SASA calculation."

    # Rename residue for bookkeeping.
    $mse_all set resname MET

    # Rename selenium atom to methionine sulfur atom and explicitly
    # assign a sulfur-like vdW radius.
    if {[$mse_se num] > 0} {
        $mse_se set name SD
        $mse_se set radius 1.80
    }

    puts "MSE -> MET conversion completed in VMD memory."
    puts "Coordinates were not changed."
    puts ""
} else {
    puts "No MSE residue found in chain $protein_chain."
    puts ""
}

$mse_all delete
$mse_se delete


# ================================================================
# Define structural environments
# ================================================================

# Protein selection.
# Include MET explicitly so the converted MSE residues are retained
# even if VMD's internal 'protein' flag was established on file load.
set protein_seltext \
    "chain $protein_chain and (protein or resname MET)"

# DNA selection.
# 'nucleic' catches standard DNA/RNA residues; explicit residue names
# retain common modified/alternative DNA residue names.
set dna_seltext \
    "chain $dna_chains and (nucleic or resname DA DT DG DC DI DU A T G C I U 5CM 5MC M5C OMC)"

# Protein-DNA complex only; excludes crystal waters and unrelated ions.
set complex_seltext \
    "(($protein_seltext) or ($dna_seltext))"

set protein [atomselect $molid $protein_seltext]
set dna     [atomselect $molid $dna_seltext]
set complex [atomselect $molid $complex_seltext]

puts "Selected protein atoms : [$protein num]"
puts "Selected DNA atoms     : [$dna num]"
puts "Selected complex atoms : [$complex num]"
puts ""

if {[$protein num] == 0} {
    puts stderr "ERROR: Protein selection is empty."
    exit 1
}

if {[$dna num] == 0} {
    puts stderr "ERROR: DNA selection is empty."
    puts stderr "Check DNA chain IDs and residue names."
    exit 1
}


# ================================================================
# Obtain unique protein residues in the requested MBD range
# ================================================================
#
# VMD's internal 'residue' index is used for unique residue identity.
# The PDB 'resid' number is retained for output.
# ================================================================

set range_sel [atomselect $molid \
    "($protein_seltext) and resid $residue_start to $residue_end"]

set residue_indices [lsort -integer -unique [$range_sel get residue]]

if {[llength $residue_indices] == 0} {
    puts stderr "ERROR: No residues found in requested range."
    exit 1
}

# Build {PDB_resid internal_residue_index} pairs and sort by PDB residue.
set residue_pairs {}

foreach ridx $residue_indices {

    set rsel [atomselect $molid \
        "chain $protein_chain and residue $ridx"]

    if {[$rsel num] == 0} {
        $rsel delete
        continue
    }

    set pdb_resid [lindex [lsort -unique [$rsel get resid]] 0]

    if {[string is integer -strict $pdb_resid]} {
        lappend residue_pairs [list $pdb_resid $ridx]
    } else {
        puts "WARNING: non-integer resid '$pdb_resid' skipped."
    }

    $rsel delete
}

set residue_pairs [lsort -integer -index 0 $residue_pairs]

puts "Residues with coordinates in requested range: [llength $residue_pairs]"
puts ""


# ================================================================
# Open output CSV
# ================================================================

set fout [open $output_csv w]

puts $fout \
"Chain,PDB_residue_number,Residue_name,Residue_label,SASA_isolated_A2,SASA_protein_A2,SASA_complex_A2,Exposure_ratio,Delta_SASA_DNA_A2,DNA_buried_ratio,Structural_category"


# Counters
set n_interior 0
set n_surface 0
set n_interface 0


# Amino acid 3-letter to 1-letter mapping
array set aa1 {
    ALA A
    ARG R
    ASN N
    ASP D
    CYS C
    GLN Q
    GLU E
    GLY G
    HIS H
    ILE I
    LEU L
    LYS K
    MET M
    PHE F
    PRO P
    SER S
    THR T
    TRP W
    TYR Y
    VAL V
}


# ================================================================
# Per-residue SASA calculation
# ================================================================
#
# A) Isolated-residue SASA:
#       measure sasa probe residue_selection
#
#    Only atoms of that residue are present in the SASA environment.
#
# B) SASA of the same residue in the full protein:
#       measure sasa probe protein -restrict residue_selection
#
#    The entire protein provides occlusion, but only the target
#    residue's accessible surface is returned.
#
# C) SASA of the same residue in protein-DNA complex:
#       measure sasa probe complex -restrict residue_selection
#
#    Protein + DNA provide occlusion, while only the target
#    protein residue's accessible surface is returned.
# ================================================================

foreach pair $residue_pairs {

    set pdb_resid [lindex $pair 0]
    set ridx      [lindex $pair 1]

    set rsel [atomselect $molid \
        "chain $protein_chain and residue $ridx"]

    if {[$rsel num] == 0} {
        continue
    }

    set resname [lindex [lsort -unique [$rsel get resname]] 0]

    if {[info exists aa1($resname)]} {
        set oneletter $aa1($resname)
    } else {
        set oneletter X
    }

    set residue_label "${oneletter}${pdb_resid}"

    # ------------------------------------------------------------
    # 1. SASA of isolated residue
    # ------------------------------------------------------------
    set sasa_isolated \
        [measure sasa $probe_radius $rsel -samples $sasa_samples]

    # ------------------------------------------------------------
    # 2. SASA contribution of this residue in complete protein
    # ------------------------------------------------------------
    set sasa_protein \
        [measure sasa $probe_radius $protein \
            -restrict $rsel \
            -samples $sasa_samples]

    # ------------------------------------------------------------
    # 3. SASA contribution of this residue in protein-DNA complex
    # ------------------------------------------------------------
    set sasa_complex \
        [measure sasa $probe_radius $complex \
            -restrict $rsel \
            -samples $sasa_samples]

    # ------------------------------------------------------------
    # Exposure ratio
    # ------------------------------------------------------------
    if {$sasa_isolated > 0.0} {
        set exposure_ratio \
            [expr {$sasa_protein / $sasa_isolated}]
    } else {
        set exposure_ratio 0.0
    }

    # ------------------------------------------------------------
    # DNA-induced burial
    #
    # Delta_SASA_DNA =
    #     SASA_protein - SASA_complex
    #
    # DNA_buried_ratio =
    #     (SASA_protein - SASA_complex) / SASA_protein
    #
    # The FINAL interface classification uses DNA_buried_ratio,
    # not the absolute Delta SASA.
    # ------------------------------------------------------------
    set delta_sasa \
        [expr {$sasa_protein - $sasa_complex}]

    if {$sasa_protein > 0.0} {
        set dna_buried_ratio \
            [expr {$delta_sasa / $sasa_protein}]
    } else {
        set dna_buried_ratio 0.0
    }

    # ------------------------------------------------------------
    # Final mutually exclusive structural classification
    # ------------------------------------------------------------
    set category \
        [classify_category \
            $dna_buried_ratio \
            $exposure_ratio \
            $dna_buried_ratio_cutoff \
            $exposure_cutoff]

    if {$category eq "Protein interior"} {
        incr n_interior
    } elseif {$category eq "Protein surface"} {
        incr n_surface
    } elseif {$category eq "Protein-DNA interface"} {
        incr n_interface
    }

    # Quote category because it contains spaces.
    puts $fout [format \
        "%s,%s,%s,%s,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,\"%s\"" \
        $protein_chain \
        $pdb_resid \
        $resname \
        $residue_label \
        $sasa_isolated \
        $sasa_protein \
        $sasa_complex \
        $exposure_ratio \
        $delta_sasa \
        $dna_buried_ratio \
        $category]

    puts [format \
        "%-5s  iso=%9.3f  protein=%9.3f  complex=%9.3f  ER=%6.3f  dSASA=%8.3f  DNA_BR=%7.4f  %s" \
        $residue_label \
        $sasa_isolated \
        $sasa_protein \
        $sasa_complex \
        $exposure_ratio \
        $delta_sasa \
        $dna_buried_ratio \
        $category]

    $rsel delete
}

close $fout


# ================================================================
# Write summary
# ================================================================

set ntotal [expr {$n_interior + $n_surface + $n_interface}]

if {$ntotal > 0} {
    set p_interior  [expr {100.0 * $n_interior  / $ntotal}]
    set p_surface   [expr {100.0 * $n_surface   / $ntotal}]
    set p_interface [expr {100.0 * $n_interface / $ntotal}]
} else {
    set p_interior 0.0
    set p_surface 0.0
    set p_interface 0.0
}

set fsum [open $output_summary w]

puts $fsum "3C2I MeCP2 MBD SASA structural classification using VMD"
puts $fsum "============================================================"
puts $fsum "Protein chain: $protein_chain"
puts $fsum "DNA chains: $dna_chains"
puts $fsum "Requested residue range: $residue_start-$residue_end"
puts $fsum "Probe radius: $probe_radius A"
puts $fsum "SASA samples: $sasa_samples"
puts $fsum ""
puts $fsum "Exposure_ratio = SASA_protein / SASA_isolated"
puts $fsum "Delta_SASA_DNA = SASA_protein - SASA_complex"
puts $fsum "DNA_buried_ratio = (SASA_protein - SASA_complex) / SASA_protein"
puts $fsum ""
puts $fsum "Classification:"
puts $fsum "DNA_buried_ratio > $dna_buried_ratio_cutoff -> Protein-DNA interface"
puts $fsum "Otherwise Exposure_ratio < $exposure_cutoff -> Protein interior"
puts $fsum "Otherwise -> Protein surface"
puts $fsum "Priority: Protein-DNA interface > Protein interior > Protein surface"
puts $fsum ""
puts $fsum [format "Protein interior      : %d (%.2f%%)" $n_interior $p_interior]
puts $fsum [format "Protein surface       : %d (%.2f%%)" $n_surface $p_surface]
puts $fsum [format "Protein-DNA interface : %d (%.2f%%)" $n_interface $p_interface]
puts $fsum [format "Total                  : %d" $ntotal]
puts $fsum ""
puts $fsum "MSE treatment:"
puts $fsum "MSE was treated as MET in VMD memory; SE was renamed SD and radius set to 1.80 A."
puts $fsum "The original PDB coordinates and the PDB file on disk were not modified."

close $fsum


# ================================================================
# Cleanup and terminal summary
# ================================================================

$range_sel delete
$protein delete
$dna delete
$complex delete

puts ""
puts "============================================================"
puts "Classification completed"
puts "============================================================"
puts "Final rules:"
puts "  DNA_buried_ratio > $dna_buried_ratio_cutoff -> Protein-DNA interface"
puts "  Otherwise Exposure_ratio < $exposure_cutoff -> Protein interior"
puts "  Otherwise -> Protein surface"
puts ""
puts [format "Protein interior      : %d (%.2f%%)" $n_interior $p_interior]
puts [format "Protein surface       : %d (%.2f%%)" $n_surface $p_surface]
puts [format "Protein-DNA interface : %d (%.2f%%)" $n_interface $p_interface]
puts [format "Total                  : %d" $ntotal]
puts ""
puts "Output:"
puts "  $output_csv"
puts "  $output_summary"
puts ""
puts "Done."

# For command-line text mode, uncomment the next line if you want
# VMD to exit automatically after finishing:
# quit
