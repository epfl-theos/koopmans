"""Input parameters for kcpw.x

For the moment, only a subset of the keywords are uncommented and can therefore be set by the user.
The majority are untested and commented out for safety."""

from pathlib import Path

from pydantic import Field
from typing_extensions import Literal

from koopmans.pydantic import BaseModel


class ControlNamelist(BaseModel):
    calculation: Literal["scf", "nscf", "bands", "relax", "md", "vc-relax",
                         "vc-md"] = Field("scf", description="A string describing the task to be performed.")
    title: str | None = Field(None, description="reprinted on output.")
    verbosity: Literal["high", "low"] = Field(
        "low",
        description="Currently two verbosity levels are implemented:  'debug' and 'medium' have the same effect as 'high'; 'default' and 'minimal' as 'low",
    )
    restart_mode: Literal["from_scratch", "restart"] = "restart"
    nstep: int = 50
    iprint: int = 10
    isave: int = 100
    tstress: bool = False
    tprnfor: bool = False
    # tabps: bool = False
    dt: float = 1.0
    ndr: int = 50
    ndw: int = 50
    outdir: Path = Field(default=".")
    prefix: str = "cp"
    pseudo_dir: Path = Field(default=".")
    # refg: float = 0.05
    max_seconds: float = 1.E+7
    ekin_conv_thr: float = 1.E-6
    etot_conv_thr: float = 1.E-4
    forc_conv_thr: float = 1.E-3
    disk_io: Literal['default'] = 'default'
    evc_restart: bool = False
    # dipfield: bool = False
    # lberry: bool = False
    # gdir: int = 0
    # nppstr: int = 0
    wf_collect: bool = False
    # lelfield: bool = False
    # nberrycyc: int = 1
    # lkpoint_dir: bool = True
    saverho: bool = True
    write_hr: bool = False
    print_real_space_density: bool = False


class SystemNamelist(BaseModel):
    """Valid keywords for the &SYSTEM namelist in kcp.x."""
    ibrav: int = -1
    celldm: dict[int, float] = Field(default_factory=dict)
    # a: float = 0.0
    # b: float = 0.0
    # c: float = 0.0
    # cosab: float = 0.0
    # cosac: float = 0.0
    # cosbc: float = 0.0
    nat: int = 0
    ntyp: int = 0
    nbnd: int = 0
    nelec: float = 0.0
    tot_charge: float = 0.0
    tot_magnetization: float | None = None
    # multiplicity: int = 0
    ecutwfc: float = 0.0
    ecutrho: float = 0.0
    # Nr1: int = 0
    # Nr2: int = 0
    # Nr3: int = 0
    # Nr1s: int = 0
    # Nr2s: int = 0
    # Nr3s: int = 0
    # nr1b: int | None = None
    # nr2b: int | None = None
    # nr3b: int | None = None
    occupations: Literal['fixed'] = 'fixed'
    # smearing:  Literal['gaussian'] = 'gaussian'
    # degauss: float = 0.0
    nelup: float = 0.0
    neldw: float = 0.0
    nspin: int = 1
    nosym: bool = False
    # nosym_evc: bool = False
    # force_symmorphic: bool = False
    # noinv: bool = False
    # ecfixed: float = 0.0
    # qcutz: float = 0.0
    # q2sigma: float = 0.01
    # input_dft: Literal['none'] = 'none'
    starting_magnetization: list[float] = Field(default_factory=list)
    # lda_plus_u: bool = False
    # hubbard_u: float = 0.0
    # edir: int = 1
    # emaxpos: float = 0.5
    # eopreg: float = 0.1
    # eamp: float = 0.0
    # la2F: bool = False
    # lspinorb: bool = False
    # noncolin: bool = False
    # lambda_: float = 1.0
    # constrained_magnetization: Literal['none'] = 'none'
    # fixed_magnetization: float = 0.0
    # B_field: float = 0.0
    assume_isolated: Literal['none'] = 'none'
    # spline_ps: bool = False
    # real_space: bool = False
    # london: bool = False
    # london_s6: float = 0.75
    # london_rcut: float = 200.0
    # do_efield: bool = False
    # ampfield: float = 0.0
    # draw_pot: bool = False
    # sortwfc_spread: bool = False
    # pot_number: int = 1
    odd_nkscalfact: bool = False
    odd_nkscalfact_empty: bool = False
    restart_odd_nkscalfact: bool = False
    wo_odd_in_empty_run: bool = False
    aux_empty_nbnd: int = 0
    restart_from_wannier_cp: bool = False
    which_file_wannier: str = " "
    wannier_empty_only: bool = False
    print_evc0_occ_empty: bool = False
    print_wfc_anion: bool = False
    index_empty_to_save: int = 1
    do_orbdep: bool = False
    do_wf_cmplx: bool = False
    do_ee: bool = False
    do_spinsym: bool = False
    f_cutoff: float = 0.01
    fixed_state: bool = False
    fixed_band: int = 1
    restart_from_wannier_pwscf: bool = False
    # impose_bloch_symm: bool = False
    # read_centers: bool = False
    # mp1: int = 1
    # mp2: int = 1
    # mp3: int = 1
    # offset_centers_occ: bool = False
    # offset_centers_emp: bool = False


class ElectronsNamelist(BaseModel):
    # emass: float = 400.0
    # emass_cutoff: float = 2.5
    # orthogonalization: Literal['ortho'] = 'ortho'
    # ortho_eps: float = 1.E-8
    # ortho_max: int = 20
    # ortho_para: int = 0
    electron_maxstep: int = 100
    # electron_dynamics: Literal['none', 'sd', 'cg', 'damp', 'verlet', 'diis'] = 'none'
    # electron_damping: float = 0.1
    # electron_velocities:  Literal['default', 'zero'] = 'default'
    # electron_temperature: Literal['not_controlled'] = 'not_controlled'
    # ekincw: float = 0.001
    # fnosee: float = 1.0
    # ampre: float = 0.0
    # grease: float = 1.0
    # startingwfc: str = 'random'
    # startingpot: str = ' '
    conv_thr: float = 1.E-6
    empty_states_maxstep: int = 100
    empty_states_ethr: float = 0.0
    # diis_size: int = 4
    # diis_nreset: int = 3
    # diis_hcut: float = 1.0
    # diis_wthr: float = 0.0
    # diis_delt: float = 0.0
    # diis_maxstep: int = 100
    # diis_rot: bool = False
    # diis_fthr: float = 0.0
    # diis_temp: float = 0.0
    # diis_achmix: float = 0.0
    # diis_g0chmix: float = 0.0
    # diis_g1chmix: float = 0.0
    # diis_nchmix: int = 3
    # diis_nrot: int = 3
    # diis_rothr: float = 0.0
    # diis_ethr: float = 0.0
    # diis_chguess: bool = False
    # mixing_mode: Literal['plain'] = 'plain'
    # mixing_fixed_ns: int = 0
    mixing_beta: float = 0.7
    # mixing_ndim: int = 8
    # diagonalization: Literal['david'] = 'david'
    # diago_thr_init: float = 0.0
    # diago_cg_maxiter: int = 20
    # diago_david_ndim: int = 4
    # diago_diis_ndim: int = 3
    # diago_full_acc: bool = False
    # sic: Literal['none'] = 'none'
    # sic_epsilon: float = 0.0
    # sic_alpha: float = 0.0
    # force_pairing: bool = False
    # fermi_energy: float = 0.0
    # n_inner: int = 2
    # niter_cold_restart: int = 1
    # lambda_cold: float = 0.03
    # rotation_dynamics: Literal['line-minimization'] = 'line-minimization'
    # occupation_dynamics: Literal['line-minimization'] = 'line-minimization'
    # rotmass: float = 0.0
    # occmass: float = 0.0
    # rotation_damping: float = 0.0
    # occupation_damping: float = 0.0
    tcg: bool = False
    maxiter: int = 100
    # passop: float = 0.3
    # niter_cg_restart: int = 20
    # etresh: float = 1.E-6
    # epol: int = 3
    # efield: float = 0.0
    # epol2: int = 3
    # efield2: float = 0.0
    # efield_cart: dict[int, float] = Field(default_factory=lambda: {1: 0.0, 2: 0.0, 3: 0.0})
    # occupation_constraints: bool = False
    do_outerloop: bool = True
    do_outerloop_empty: bool = True
    # reortho: bool = False


class IonsNamelist(BaseModel):
    # phase_space: Literal['full', 'coarse-grained'] = 'full'
    # ion_dynamics: Literal['sd', 'cg', 'damp', 'verlet', 'none', 'bfgs', 'beeman'] = 'none'
    # ion_radius: float = 0.5
    # ion_damping: float = 0.1
    # ion_positions: Literal['default', 'from_input'] = 'default'
    # ion_velocities: Literal['zero', 'default', 'from_input'] = 'default'
    # ion_temperature: Literal['not_controlled'] = 'not_controlled'
    # tempw: float = 300.0
    # fnosep: float = -1.0
    # nhpcl: int = 0
    # nhptyp: int = 0
    # ndega: int = 0
    # tranp: bool = False
    # amprp: float = 0.0
    # greasp: float = 1.0
    # tolp: float = 100.0
    # ion_nstepe: int = 1
    # ion_maxstep: int = 100
    # delta_t: float = 1.0
    # nraise: int = 1
    # refold_pos: bool = False
    # remove_rigid_rot: bool = False
    # upscale: float = 10.0
    # pot_extrapolation: Literal['atomic'] = 'atomic'
    # wfc_extrapolation: Literal['none'] = 'none'
    # num_of_images: int = 0
    # first_last_opt: bool = False
    # use_masses: bool = False
    # use_freezing: bool = False
    # opt_scheme: Literal['quick-min'] = 'quick-min'
    # temp_req: float = 0.0
    # ds: float = 1.0
    # path_thr: float = 0.05
    # ci_scheme: Literal['no-CI'] = 'no-CI'
    # k_max: float = 0.1
    # k_min: float = 0.1
    # fixed_tan: bool = False
    # bfgs_ndim: int = 1
    # trust_radius_max: float = 0.8
    # trust_radius_min: float = 1.E-4
    # trust_radius_ini: float = 0.5
    # w_1: float = 0.01
    # w_2: float = 0.50
    # sic_rloc: float = 0.0
    # fe_step: float = 0.4
    # fe_nstep: int = 100
    # sw_nstep: int = 10
    # eq_nstep: int = 0
    # g_amplitude: float = 0.005
    pass


class CellNamelist(BaseModel):
    # cell_parameters: Literal['default'] = 'default'
    # cell_dynamics: Literal['sd', 'pr', 'none', 'w', 'damp-pr', 'damp-w', 'bfgs'] = 'none'
    # cell_velocities: Literal['default', 'zero'] = 'default'
    # press: float = 0.0
    # wmass: float = 0.0
    # cell_temperature: Literal['nose', 'not_controlled', 'rescaling'] = 'not_controlled'
    # temph: float = 0.0
    # fnoseh: float = 1.0
    # greash: float = 1.0
    # cell_dofree: Literal['all', 'volume', 'x', 'y', 'z', 'xy', 'xz', 'yz', 'xyz'] = 'all'
    # cell_factor: float = 0.0
    # cell_nstepe: int = 1
    # cell_damping: float = 0.0
    # press_conv_thr: float = 0.5
    pass


class EENamelist(BaseModel):
    which_compensation: Literal['none'] = 'none'
    tcc_odd: bool = False


class NKSICNamelist(BaseModel):
    esic_conv_thr: float = 1.E-5
    # do_nk: bool = False
    do_pz: bool = False
    do_nki: bool = False
    # do_nkpz: bool = False
    do_nkipz: bool = False
    do_innerloop: bool = Field(default=False, description="main switch of inner loop minimization")
    do_innerloop_empty: bool = Field(
        default=False, description="main switch of inner loop minimization for empty states")
    # l_comp_cmplxfctn_index: bool = Field(default=False, description="compute the complexification index")
    do_innerloop_cg: bool = Field(default=False, description="main switch of cg inner loop minimization")
    innerloop_dd_nstep: int = Field(
        default=50, description="number of outer loop damped dynamics steps between each inner loop minimization")
    innerloop_cg_nsd: int = Field(
        default=20, description="number of initial steepest-descent steps in cg inner loop minimization")
    innerloop_cg_nreset: int = Field(
        default=10, description="number of cg steps after which the search direction is set to the steepest-descent direction in inner loop minimization")
    innerloop_nmax: int = Field(default=10000, description="maximum number of inner loop steps")
    innerloop_cg_ratio: float = Field(default=1.e-3)
    # innerloop_init_n: int | None = None
    # innerloop_until: int = -1
    # innerloop_atleast: int = 0
    nkscalfact: float = Field(default=1.0, description="NK coeffcient")
    # hfscalfact: float = Field(default=1.0, description="HF coefficient")
    # nknmax: int = Field(default=-1, description="if <> -1, index of the last orbital on which NK is applied")
    # do_hf: bool = Field(default=False, description="main switch for HF calculations")
    # do_wxd: bool = Field(default=True, description="include cross-terms in NK potential")
    # do_wref: bool = Field(default=True, description="include reference variational terms")
    # do_pz_renorm: bool = Field(default=False)
    do_bare_eigs: bool = Field(default=False)
    # kfact: float = 0.0
    # fref: float = 0.5
    # rhobarfact: float = 1.0
    # vanishing_rho_w: float = 1.0e-12
    which_orbdep: str | None = None
    # iprint_spreads: int = -1
    # iprint_manifold_overlap: int = -1
    hartree_only_sic: bool = False
    # finite_field_introduced: bool = False
    # finite_field_for_empty_state: bool = False


class KCPInputParameters(BaseModel):
    control: ControlNamelist = Field(default_factory=lambda: ControlNamelist())
    system: SystemNamelist = Field(default_factory=lambda: SystemNamelist())
    electrons: ElectronsNamelist = Field(default_factory=lambda: ElectronsNamelist())
    ions: IonsNamelist = Field(default_factory=lambda: IonsNamelist())
    cell: CellNamelist = Field(default_factory=lambda: CellNamelist())
    ee: EENamelist = Field(default_factory=lambda: EENamelist())
    nksic: NKSICNamelist = Field(default_factory=lambda: NKSICNamelist())
