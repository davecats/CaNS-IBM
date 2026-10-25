#!/usr/bin/env python3
# MacDonald et al. (JFM 2019) set A full-span case (Re_tau = 395, k = h/18, lambda = 7.07k) for CaNS-IBM.
# usage: make_input.py <rough|smooth> <phase: cfr|cpg> [key=value ...]   -> input.nml on stdout
# cfr: spin-up at constant bulk velocity (velf); cpg: constant pressure gradient bforce = 1 (u_tau = 1, h = 1).
import sys, math
case, phase = sys.argv[1], sys.argv[2]
o = dict(nx=512, ny=256, nlx=16, nly=8, nz=480 if case=='rough' else 240, ndims=4,
         restart='F', tmax=1000., twmax=46., iout0d=100, iout1d=1000, iout2d=0, iout3d=0, isave=5000,
         iwall=0, cfl=0.95, nstep=100000000)
for a in sys.argv[3:]:
    key, v = a.split('='); o[key] = type(o[key])(v) if key in o and not isinstance(o[key], str) else v
k = 1./18.; lam = 2*math.pi/16.
lx, ly = o['nlx']*lam, o['nly']*lam
rough = case == 'rough'
lz = 2. + 2.*k if rough else 2.
velf = 10.5 if rough else 17.4
cfr = phase == 'cfr'
T, F = 'T', 'F'
print(f"""&dns
ng(1:3) = {o['nx']}, {o['ny']}, {o['nz']}
l(1:3) = {lx!r}, {ly!r}, {lz!r}
gtype = 5, gr = 0.
dzp_fine = 0.35, dzp_coarse = 6.9
cfl = {o['cfl']}, dtmax = 1.e9, dt_f = -1.
visci = 395.
inivel = 'log'
is_wallturb = T
nstep = {o['nstep']}, time_max = {o['tmax']}, tw_max = {o['twmax']}
stop_type(1:3) = T, T, T
restart = {o['restart']}, is_overwrite_save = T, nsaves_max = 0
icheck = 100, iout0d = {o['iout0d']}, iout1d = {o['iout1d']}, iout2d = {o['iout2d']}, iout3d = {o['iout3d']}, isave = {o['isave']}
cbcvel(0:1,1:3,1) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,2) = 'P','P',  'P','P',  'D','D'
cbcvel(0:1,1:3,3) = 'P','P',  'P','P',  'D','D'
cbcpre(0:1,1:3)   = 'P','P',  'P','P',  'N','N'
bcvel(0:1,1:3,1) =  0.,0.,   0.,0.,   0.,0.
bcvel(0:1,1:3,2) =  0.,0.,   0.,0.,   0.,0.
bcvel(0:1,1:3,3) =  0.,0.,   0.,0.,   0.,0.
bcpre(0:1,1:3)   =  0.,0.,   0.,0.,   0.,0.
bforce(1:3) = {0. if cfr else 1.}, 0., 0.
is_forced(1:3) = {T if cfr else F}, F, F
velf(1:3) = {velf}, 0., 0.
nscal = 1
gacc = 0.0,0.0,0.0
dims(1:2) = {o['ndims']}, 1, ipencil_axis = 1
/

&cudecomp
cudecomp_t_comm_backend = 0, cudecomp_is_t_enable_nccl = T, cudecomp_is_t_enable_nvshmem = T
cudecomp_h_comm_backend = 0, cudecomp_is_h_enable_nccl = T, cudecomp_is_h_enable_nvshmem = T
/

&scalar
alphai(1) = 276.5
beta = 0
iniscal(1) = 'zer'
cbcscal(0:1,1:3,1) = 'P','P',  'P','P',  'D','D'
bcscal(0:1,1:3,1)  =  0.,0.,   0.,0.,   0.,0.
cbcscal_ibm(1) = 'D'
bcscal_ibm(1)  = 0.
const_scal = F
scal_mean_grad(1) = -1.
ssource(1) = 0.
is_sforced(1) = F
scalf(1) = 0.
is_boussinesq_buoyancy = F
/

&numerics
is_impdiff = T, is_impdiff_1d = T
is_poisson_pcr_tdma = F
is_fast_mom_kernels = F
adv_form = 'vv'
/

&other_options
is_debug = T, is_timing = T
/

&io
io_backend = 'mpiio'
/

&ibm
is_ibm = {T if rough else F}
ibm_2nd = {T if rough else F}
do_richardson = T
ibm_direction = F,F, F,F, T,T
l_0 = 0.0,0.0, 0.0,0.0, 0.0,{lz!r}
n_wave = {o['nlx']},{o['nlx']}, {o['nly']},{o['nly']}, 0,0
amp_l = {k!r},{k!r}, 0.0,0.0, 0.0,0.0
phase_l = 0.0,0.0, 0.0,0.0, 0.0,0.0
ibm_shape = 'macdonald'
ibm_hch = 1.
bulk_vel_def = 'intrinsic'
ibm_proj = '{'damped' if rough else 'none'}'
ibm_iout_wall = {o['iwall']}
/

&hmap
use_hmap = F
hmap_mode = "fit"
override_grid = F
hmap_name = "thakkar"
/""")
