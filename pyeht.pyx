import numpy as np
#cimport numpy as cnp

from libc.stdlib cimport malloc, calloc, free
from libc.stdio cimport FILE
from libc.string cimport strcpy

ctypedef char BOOLEAN
ctypedef double real


cdef extern from "bind.h":
    int ATOM_SYMB_LEN=10
    int FAT=6
    int THIN=1
    int MOLECULAR=27
    real THE_CONST=1.75
    real NN_DISTANCE=2.5
    real MULLER_MIX_DEF=0.20
    real MULLER_E_TOL_DEF=0.01
    real MULLER_Z_TOL_DEF=0.001
    real K_OFFSET=0.01
    int num_orbs
    int* orbital_lookup_table
    ctypedef struct Z_mat_type:
        pass
    ctypedef struct point_type:
        real x, y, z
    ctypedef struct atom_type:
        char symb[10]
        char chg_it_vary
        int which_atom, at_number, num_valence
        point_type loc
        Z_mat_type Zmat_loc
        int ns, np, nd, nf
        real coul_s, coul_p, coul_d, coul_f
        real coeff_d1, coeff_d2, coeff_f1, coeff_f2
        real s_A, s_B, s_C
        real p_A, p_B, p_C
        real d_A, d_B, d_C
        real muller_s_E[7], muller_p_E[7], muller_d_E[7]
        real muller_s_Z[4], muller_p_Z[4], muller_d_Z[4]
        real init_s_occup, init_p_occup, init_d_occup
    ctypedef struct geom_frag_type:
        pass
    ctypedef struct Tvect_type:
        pass
    ctypedef struct xtal_defn_type:
        pass
    ctypedef struct equiv_atom_type:
        pass
    ctypedef struct cell_type:
        int dim, num_atoms, num_raw_atoms
        atom_type* atoms
        geom_frag_type* geom_frags
        char using_Zmat, using_xtal_coords
        real* distance_mat
        real num_electrons, charge
        char* sym_elems
        Tvect_type tvects[3]
        point_type recip_vects[3]
        xtal_defn_type xtal_defn
        int overlaps[3]
        point_type COM
        real princ_axes[3][3]
        equiv_atom_type* equiv_atoms

    ctypedef struct printing_info_type:
        pass
    ctypedef struct chg_it_parm_type:
        pass
    ctypedef struct orbital_occup_type:
        pass
    ctypedef struct k_point_type:
        point_type loc
        real weight
        real num_filled_bands
    ctypedef struct FMO_frag_type:
        pass
    ctypedef struct FMO_prop_type:
        pass
    ctypedef struct p_DOS_type:
        pass
    ctypedef struct COOP_type:
        pass
    ctypedef struct walsh_details_type:
        int num_steps
        int num_vars
        real* values
        printing_info_type* things_to_print
    ctypedef struct band_info_type:
        pass
    ctypedef struct overlap_cancel_type:
        pass
    ctypedef struct detail_type:
        char title[240]
        char filename[240]
        int Execution_Mode
        BOOLEAN just_geom, avg_props, gradients, save_energies
        BOOLEAN use_symmetry, find_princ_axes
        BOOLEAN vary_zeta
        BOOLEAN eval_electrostat
        BOOLEAN weighted_Hij
        BOOLEAN dump_overlap, dump_hamil
        BOOLEAN dump_sparse_mats
        BOOLEAN dump_dist_mat
        # printing options
        int upper_level_PRT, lower_level_PRT
        real max_dist_PRT
        BOOLEAN distance_mat_PRT
        BOOLEAN chg_mat_PRT, Rchg_mat_PRT, wave_fn_PRT
        BOOLEAN net_chg_PRT, overlap_mat_PRT, electrostat_PRT, fermi_PRT
        BOOLEAN hamil_PRT, energies_PRT, levels_PRT
        BOOLEAN avg_OP_mat_PRT, avg_ROP_mat_PRT
        BOOLEAN no_total_DOS_PRT
        BOOLEAN just_avgE, just_matrices
        BOOLEAN mod_OP_mat_PRT, mod_ROP_mat_PRT, mod_net_chg_PRT
        BOOLEAN orbital_mapping_PRT
        int line_width
        char diag_wo_overlap
        char store_R_overlaps
        printing_info_type* step_print_options
        int num_MOs_to_print
        int *MOs_to_print
        real rho
        char do_chg_it
        chg_it_parm_type chg_it_parms
        real close_nn_contact
        int num_occup_AVG
        real occup_AVG_step
        int num_orbital_occups
        orbital_occup_type* orbital_occups
        int num_KPOINTS
        k_point_type* K_POINTS
        char use_automatic_kpoints, use_high_symm_p
        int points_per_axis[3]
        real k_offset
        int num_occup_KPOINTS
        real* occup_KPOINTS
        int num_bonds_OOP
        int num_FMO_frags
        int num_FCO_frags
        FMO_frag_type *FMO_frags
        FMO_prop_type *FMO_props
        int num_proj_DOS
        p_DOS_type *proj_DOS
        COOP_type *the_COOPS
        char do_moments
        int num_moments
        real* moments
        walsh_details_type walsh_details
        band_info_type *band_info
        int num_overlaps_off
        overlap_cancel_type *overlaps_off
        real the_const
        real sparsify_value
        real symm_tol
        int num_sym_ops
        real *characters
        char do_muller_it
        int *atoms_to_vary
        real muller_mix, muller_E_tol, muller_Z_tol

    void charge_to_num_electrons(cell_type*)

cdef extern from "symmetry.h":
    real SYMM_TOL=1e-3

cdef extern from "prototypes.h":
    void fill_atomic_parms(atom_type*, int, FILE*, char*)
    void build_orbital_lookup_table(cell_type*, int*, int**)
    void check_for_errors(cell_type*, detail_type*, int)

"""
need to write function which replicates functionality of 'run_bind'

- allocate cell_type
- allocate detail_type
- set FILE handles (ie dev/null/)?
- read_inputfile
- automagic_k_points???
- check_for_errors
- inner_wrapper
- steal results!
  - reassign pointer to results memory to something else
  - ie avoid copying memory?
- cleanup_memory
- free cell_type and detail_type

"""
cdef detail_type* make_details():
    """Make a default version of a details object"""
    cdef detail_type* details

    details = <detail_type*> malloc(sizeof(detail_type))

    # Initial values
    details.walsh_details.num_steps = 1
    details.walsh_details.num_vars = 0
    details.use_symmetry = 0
    details.find_princ_axes = 0
    details.vary_zeta = 0
    details.avg_props = 0
    details.just_geom = 0
    details.dump_overlap = 0
    details.dump_hamil = 0
    details.sparsify_value = 0.0
    details.Execution_Mode = FAT
    details.the_const = THE_CONST
    details.weighted_Hij = 1
    details.eval_electrostat = 0
    details.close_nn_contact = NN_DISTANCE
    details.symm_tol = SYMM_TOL
    details.muller_mix = MULLER_MIX_DEF
    details.muller_E_tol = MULLER_E_TOL_DEF
    details.muller_Z_tol = MULLER_Z_TOL_DEF
    details.num_moments = 4
    details.line_width = 80
    details.k_offset = K_OFFSET

    return details

cdef void customise_details(detail_type* details):
    """settings for how we run a tight binding calculation"""
    # MOLECULAR run type
    details.Execution_Mode = MOLECULAR
    details.num_KPOINTS = 1
    details.K_POINTS = <k_point_type*>calloc(1,sizeof(k_point_type))
    details.K_POINTS[0].weight = 1.0
    # charge 0 -> cell

    # Nonweighted
    details.weighted_Hij = 0

    # dump hamiltonian
    details.dump_hamil = 1
    # dump overlap
    details.dump_overlap = 1
    # Just Matrices
    details.just_matrices = 1


cdef cell_type* make_cell():
    cdef cell_type* cell

    cell = <cell_type*> malloc(sizeof(cell_type))
    # initial values
    cell.equiv_atoms = NULL  # isn't an int, but gets set to 0?
    cell.charge = -1000.0

    return cell

cdef void customise_cell(cell_type* cell, int num_atoms,
                         double[:, ::1] positions,
                         char[:, ::1] atom_types,
                         double charge):
    """Set the cell to our simulation contents

    Parameters
    ----------
    cell : ptr to cell_type
      the cell object to modify in place
    num_atoms : int
      size of simulation
    positions : double
      positions of all atoms
    atom_type : char
      names of atoms, encoded into int8 array
    charge : float
      net charge of system
    """
    cdef int i
    cdef bytes name

    # load geometry into program
    cell.using_Zmat = 0
    cell.using_xtal_coords = 0

    cell.num_atoms = num_atoms
    cell.atoms = <atom_type*>calloc(cell.num_atoms, sizeof(atom_type))
    for i in range(num_atoms):
        strcpy(cell.atoms[i].symb, &atom_types[i][0])
        cell.atoms[i].loc.x = positions[i][0]
        cell.atoms[i].loc.y = positions[i][1]
        cell.atoms[i].loc.z = positions[i][2]
        cell.atoms[i].which_atom = i
    cell.num_raw_atoms = cell.num_atoms

    # Figure out parameters
    fill_atomic_parms(cell.atoms, cell.num_atoms, NULL, NULL)

    # Charge 0
    cell.geom_frags = NULL
    cell.num_electrons = 0
    cell.charge = 0
    charge_to_num_electrons(cell)


def run_eht(double[:, ::1] positions, list elements, double charge):
    """Run tight binding calculations

    Parameters
    ----------
    positions
    elements

    Returns
    -------
    Hamiltonian
    Overlap
    """
    cdef cell_type* cell
    cdef detail_type* details
    cdef int num_atoms
    cdef char[:, ::1] atom_types

    # TODO Sort out file handles here

    num_atoms = positions.shape[0]

    # TODO Add check for null terminated strings
    atom_types = np.array(elements, dtype='S10').view(
        np.int8).reshape(num_atoms, -1)

    cell = make_cell()
    details = make_details()
    customise_details(details)
    customise_cell(cell, num_atoms, positions, atom_types, charge)

    #build_orbital_lookup_table(cell, &num_orbs, &orbital_lookup_table)

    n = cell.num_electrons

    free(cell)
    free(details)

    return n
