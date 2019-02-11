import cython

import numpy as np
#cimport numpy as cnp

from libc.stdlib cimport malloc, calloc, free
from libc.stdio cimport FILE, fopen, stdout
from libc.string cimport strcpy
from libc.signal cimport signal, SIGINT

## cnp.import_array()
## define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

ctypedef char BOOLEAN
ctypedef double real_t


cdef extern from "bind.h":
    FILE *status_file, *output_file, *walsh_file, *band_file, *FMO_file
    FILE *MO_file
    int ATOM_SYMB_LEN=10
    int FAT=6
    int THIN=1
    int MOLECULAR=27
    real_t THE_CONST=1.75
    real_t NN_DISTANCE=2.5
    real_t MULLER_MIX_DEF=0.20
    real_t MULLER_E_TOL_DEF=0.01
    real_t MULLER_Z_TOL_DEF=0.001
    real_t K_OFFSET=0.01
    cell_type* unit_cell
    detail_type* details
    int num_orbs
    int* orbital_lookup_table
    ctypedef struct hermetian_matrix_type:
        int dim
        real_t *mat
    hermetian_matrix_type Hamil_R, Hamil_K, Overlap_R, Overlap_K
    ctypedef struct Z_mat_type:
        pass
    ctypedef struct point_type:
        real_t x, y, z
    ctypedef struct atom_type:
        char symb[10]
        char chg_it_vary
        int which_atom, at_number, num_valence
        point_type loc
        Z_mat_type Zmat_loc
        int ns, np, nd, nf
        real_t coul_s, coul_p, coul_d, coul_f
        real_t coeff_d1, coeff_d2, coeff_f1, coeff_f2
        real_t s_A, s_B, s_C
        real_t p_A, p_B, p_C
        real_t d_A, d_B, d_C
        real_t muller_s_E[7], muller_p_E[7], muller_d_E[7]
        real_t muller_s_Z[4], muller_p_Z[4], muller_d_Z[4]
        real_t init_s_occup, init_p_occup, init_d_occup
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
        real_t* distance_mat
        real_t num_electrons, charge
        char* sym_elems
        Tvect_type tvects[3]
        point_type recip_vects[3]
        xtal_defn_type xtal_defn
        int overlaps[3]
        point_type COM
        real_t princ_axes[3][3]
        equiv_atom_type* equiv_atoms
    ctypedef struct printing_info_type:
        pass
    ctypedef struct chg_it_parm_type:
        pass
    ctypedef struct orbital_occup_type:
        pass
    ctypedef struct k_point_type:
        point_type loc
        real_t weight
        real_t num_filled_bands
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
        real_t* values
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
        real_t max_dist_PRT
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
        real_t rho
        char do_chg_it
        chg_it_parm_type chg_it_parms
        real_t close_nn_contact
        int num_occup_AVG
        real_t occup_AVG_step
        int num_orbital_occups
        orbital_occup_type* orbital_occups
        int num_KPOINTS
        k_point_type* K_POINTS
        char use_automatic_kpoints, use_high_symm_p
        int points_per_axis[3]
        real_t k_offset
        int num_occup_KPOINTS
        real_t* occup_KPOINTS
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
        real_t* moments
        walsh_details_type walsh_details
        band_info_type *band_info
        int num_overlaps_off
        overlap_cancel_type *overlaps_off
        real_t the_const
        real_t sparsify_value
        real_t symm_tol
        int num_sym_ops
        real_t *characters
        char do_muller_it
        int *atoms_to_vary
        real_t muller_mix, muller_E_tol, muller_Z_tol

    void charge_to_num_electrons(cell_type*)

cdef extern from "symmetry.h":
    real_t SYMM_TOL=1e-3

cdef extern from "prototypes.h":
    void fill_atomic_parms(atom_type*, int, FILE*, char*)
    void build_orbital_lookup_table(cell_type*, int*, int**)
    void check_for_errors(cell_type*, detail_type*, int)
    void run_eht(FILE*)
    void inner_wrapper(char*, bool)
    void cleanup_memory()

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
cdef void make_details(detail_type* details):
    """Make a default version of a details object"""
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

    # set other zeros
    # details.title = 'pyeht'
    details.num_FMO_frags = 0
    details.num_FCO_frags = 0
    details.do_chg_it = 0
    details.use_automatic_kpoints = 0
    details.use_high_symm_p = 0
    details.no_total_DOS_PRT = 0
    details.num_proj_DOS = 0
    details.the_COOPS = NULL
    details.num_MOs_to_print = 0
    details.band_info = NULL
    details.just_avgE = 0
    details.just_matrices = 0
    details.do_muller_it = 0


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
    #details.dump_hamil = 1
    # dump overlap
    #details.dump_overlap = 1
    # Just Matrices
    details.just_matrices = 1


cdef void make_cell(cell_type* cell):
    # initial values
    cell.equiv_atoms = NULL  # isn't an int, but gets set to 0?
    cell.geom_frags = NULL
    cell.charge = -1000.0
    cell.using_Zmat = 0
    cell.using_xtal_coords = 0


@cython.boundscheck(False)
@cython.wraparound(False)
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

cdef void steal_matrix(real_t* mat, int n, double[::1] newmat):
    """Does something a little like dump_hermetian_mat"""
    #mat = np.array(mat)
    cdef int i

    for i in range(n):
        newmat[i] = mat[i]


def run_yaehmop(double[:, ::1] positions, list elements, double charge):
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
    #cdef cell_type* my_cell
    #cdef detail_type* my_details
    cdef int num_atoms
    cdef char[:, ::1] atom_types
    #cdef int num_orbs
    #cdef int* orbital_lookup_table
    global orbital_lookup_table
    global num_orbs
    cdef FILE *hell
    global status_file
    global unit_cell
    global details
    # results
    global Hamil_R, Hamil_K, Overlap_R, Overlap_K


    details = <detail_type*> calloc(1, sizeof(detail_type))
    unit_cell = <cell_type*> calloc(1, sizeof(cell_type))

    # file handles (to the pits of hell)
    hell = fopen('/dev/null', 'w')
    #hell = stdout
    status_file = hell
    output_file = hell
    walsh_file = hell
    band_file = hell
    FMO_file = hell
    MO_file = hell

    num_atoms = positions.shape[0]
    # TODO Add check for null terminated strings
    atom_types = np.array(elements, dtype='S10').view(
        np.int8).reshape(num_atoms, -1)

    make_cell(unit_cell)
    make_details(details)

    customise_details(details)
    customise_cell(unit_cell, num_atoms, positions, atom_types, charge)
    build_orbital_lookup_table(unit_cell, &num_orbs, &orbital_lookup_table)

    run_eht(hell)
    #inner_wrapper('no-file', 1)

    H_mat = np.empty(num_orbs * num_orbs, dtype=np.float64)
    S_mat = np.empty(num_orbs * num_orbs, dtype=np.float64)

    #steal_matrix(Hamil_R.mat, num_orbs * num_orbs, H_mat)
    #steal_matrix(Overlap_R.mat, num_orbs * num_orbs, S_mat)

    # Once we're done grabbing results, free memory again
    cleanup_memory()
    # These objects are cleaned last now their contents are free
    free(unit_cell)
    free(details)

    return H_mat, S_mat
