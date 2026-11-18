#!/bin/bash
# =====================================================================================================
#  N2 / 6-31G, R = 3.0 bohr, N 1s frozen, point-group projection C4 (irrep A, from Sigma_g+)
#  frozen core -> RHF / UHF -> C4 x SUHF -> C4 x SUCCSD (2 x 100 cycles) -> coherent (T) and +(T)
# =====================================================================================================
#  Run one step at a time and check its log before going on:
#     bash run.sh 1    frozen core, integrals, C4v projector cut to C4      -> FC.p (C4), FC_C4v.p
#     bash run.sh 2    RHF / UHF in the frozen-core basis (Main_HF.py)      -> GHFMO.p (stable UHF)
#     bash run.sh 3    gradient C4 x SUHF, one basin hop, from the UHF      -> SGHF.p
#     bash run.sh 4    Fock-based C4 x SUHF from the gradient solution      -> SGHFFock.p (orbitals, eff.-Fock energies)
#     bash run.sh 5    C4 x SUCCSD, first 100 Broyden cycles                -> pcc/
#     bash run.sh 6    C4 x SUCCSD, 100 more cycles (restart in pcc/)       -> pcc/
#     bash run.sh 7    coherent (T) and +(T)                                -> pcc/result_T_coh.json, result_T_plus.json
#  The C4 x SUCCSD residual decreases slowly: 2 x 100 cycles bring it to about 5e-4 (energy settled to about
#  1e-4 Eh), not to the default 5e-5. Restarting the solver (step 6) helps; run step 6 again for more cycles.
#  Settings: parameter.py (molecule, frozen core, point group C4v / A1), parameter_fc.py (spin projection).
#  The PCC module must be compiled first (README, "Compiling the code").
# =====================================================================================================
set -eo pipefail
cd "$(dirname "$0")"

REPO=$(cd ../.. && pwd)
export PYTHONPATH=$REPO
export LD_LIBRARY_PATH=$REPO/SUCCSD_new:$LD_LIBRARY_PATH

# spin projection of the PCC steps (as in parameter_fc.py); the point-group projector is read from FC.p
SP=2
GRID="1 8"


# ----------------------------------------------------------------------------------------------- helper
# Run a script of the repository so that it imports parameter.py / parameter_fc.py of THIS folder
# ("python -c" puts the current folder first on the import path).
run_repo_script() {
    local script=$1
    shift
    python -c "import sys; sys.argv = ['$script'] + sys.argv[1:]; exec(open('$REPO/$script').read())" "$@"
}

# C4 x SUCCSD in the effective-Fock orbitals, at most 100 Broyden cycles; extra options are passed on
run_pccsd_100() {
    local log=$1
    shift
    PCC_MOE_FILE=$PWD/moe.txt PCC_CYCMAX=100 python $REPO/SUCCSD_new/run_pcc.py \
        --module PUCC_X2 --libdir $REPO/SUCCSD_new \
        --fc FC.p --mo SGHFFock.p --sp $SP --grid $GRID --pg 1 \
        --semicanon 0 --work pcc "$@" | tee $log
}


# ----------------------------------------------------------------------------------------------- steps
step1_frozen_core() {
    run_repo_script Main_FC.py | tee step1_fc.log           # C4v projector (PG = 'C4v' in parameter.py)
    cp FC.p FC_C4v.p
    # C4 = rotation subgroup of C4v: keep the first four operators (E, C4, C2, C4^3) and their characters
    python -c "
import pickle
fc = pickle.load(open('FC.p', 'rb'))
fc[4] = fc[4][:4]      # characters (weightpg)
fc[5] = fc[5][:4]      # operators (Rpg)
pickle.dump(fc, open('FC.p', 'wb'), protocol=4)
print('C4 projector: characters', fc[4][:, 0, 0].real)
" | tee -a step1_fc.log
}

step2_hf() {
    run_repo_script Main_HF.py | tee step2_hf.log
}

step3_gradient_phf() {
    rm -f SGHF.p hubsghf.p                    # restart files of an earlier PHF run
    python -c "
import PHF
PHF.nhop = 1           # one basin hop from the UHF (0: plain BFGS, stays at the nearest stationary point)
PHF.fixgauge = False   # SUHF: no spin-gauge fixing
exec(open('$REPO/Main_PHF.py').read())
" 2>&1 | tee step3_phf.log
    if [ ! -f SGHF.p ]; then
        # the Fortran code stops without an error status (e.g. "Not Hermitian" after a random basin-hop step)
        echo "STOP: the gradient PHF did not finish (see step3_phf.log); rerun this step"
        exit 1
    fi
    # the basin hop is random: check the projected norm 'Ovlp' in the log (here about 0.64); a value near zero
    # means the projection collapsed, rerun this step
}

step4_fock_phf() {
    run_repo_script PHFFock.py SGHF.p | tee step4_fock.log
    if grep -q "converged: True" step4_fock.log; then
        return
    fi
    echo "Fock PHF not converged: retry once from the lowest-energy iterate"
    run_repo_script PHFFock.py newmo_best.p | tee step4_fock_retry.log
    if ! grep -q "converged: True" step4_fock_retry.log; then
        echo "STOP: no stationary PHF solution, do not run PCC on it"
        exit 1
    fi
}

step5_pccsd_first_100() {
    # effective-Fock orbital energies (5th entry of SGHFFock.p) for the preconditioner
    python -c "import pickle, numpy; numpy.savetxt('moe.txt', pickle.load(open('SGHFFock.p', 'rb'))[4])"
    run_pccsd_100 step5_pcc.log
}

step6_pccsd_next_100() {
    # --keep-chk: start from the amplitudes of the previous run (Fortran restart files Chk_PUCC_* in pcc/)
    run_pccsd_100 step6_pcc.log --keep-chk
}

step7_triples() {
    local e_pcc
    e_pcc=$(python -c "import json; print(json.load(open('pcc/result.json'))['E_PCC'])")

    # coherent (T): projected (T) in the CC orbitals, effective-Fock orbital energies, with the C4 projector
    python $REPO/SUCCSD_new/run_pcc_t.py \
        --fc FC.p --mo pcc/MO_used.p --t1t2 pcc/T1T2.p \
        --newmo SGHFFock.p --moe eff \
        --sp $SP --grid $GRID --pg 1 --epcc $e_pcc --tag _coh | tee step7_coh_t.log

    # +(T): semicanonical orbitals of the bare Fock in the CC orbital spaces, then single-determinant (T)
    python $REPO/SUCCSD_new/semicanon_plusT.py FC.p pcc/MO_used.p $SP MO_sc.p | tee step7_semicanon.log
    python $REPO/SUCCSD_new/run_pcc_t.py \
        --fc FC.p --mo pcc/MO_used.p --t1t2 pcc/T1T2.p \
        --newmo MO_sc.p --moe fock \
        --sp 0 --grid 1 1 --epcc $e_pcc --tag _plus | tee step7_plus_t.log
}


# ----------------------------------------------------------------------------------------------- main
case "$1" in
    1) step1_frozen_core ;;
    2) step2_hf ;;
    3) step3_gradient_phf ;;
    4) step4_fock_phf ;;
    5) step5_pccsd_first_100 ;;
    6) step6_pccsd_next_100 ;;
    7) step7_triples ;;
    *) sed -n '2,19p' "$0"; exit 1 ;;
esac
