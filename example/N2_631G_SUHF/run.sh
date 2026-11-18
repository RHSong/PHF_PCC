#!/bin/bash
# =====================================================================================================
#  N2 / 6-31G, R = 3.0 bohr, N 1s frozen
#  frozen core -> RHF / UHF -> SUHF -> SUCCSD -> coherent (T) and +(T)
# =====================================================================================================
#  Run one step at a time and check its log before going on:
#     bash run.sh 1    frozen core and integrals (Main_FC.py)               -> FC.p
#     bash run.sh 2    RHF / UHF in the frozen-core basis (Main_HF.py)      -> GHFMO.p (stable UHF)
#     bash run.sh 3    gradient SUHF, one basin hop, from the UHF           -> SGHF.p
#     bash run.sh 4    Fock-based SUHF from the gradient solution           -> SGHFFock.p (orbitals, eff.-Fock energies)
#     bash run.sh 5    SUCCSD in the effective-Fock orbitals                -> pcc/
#     bash run.sh 6    coherent (T) and +(T)                                -> pcc/result_T_coh.json, result_T_plus.json
#  Settings: parameter.py (molecule, frozen core, point group), parameter_fc.py (spin projection).
#  The PCC module must be compiled first (README, "Compiling the code").
# =====================================================================================================
set -eo pipefail
cd "$(dirname "$0")"

REPO=$(cd ../.. && pwd)
export PYTHONPATH=$REPO
export LD_LIBRARY_PATH=$REPO/SUCCSD_new:$LD_LIBRARY_PATH

# spin projection of the PCC steps (as in parameter_fc.py)
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


# ----------------------------------------------------------------------------------------------- steps
step1_frozen_core() {
    run_repo_script Main_FC.py | tee step1_fc.log
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

step5_pccsd() {
    # effective-Fock orbital energies (5th entry of SGHFFock.p) for the preconditioner
    python -c "import pickle, numpy; numpy.savetxt('moe.txt', pickle.load(open('SGHFFock.p', 'rb'))[4])"
    PCC_MOE_FILE=$PWD/moe.txt python $REPO/SUCCSD_new/run_pcc.py \
        --module PUCC_X2 --libdir $REPO/SUCCSD_new \
        --fc FC.p --mo SGHFFock.p --sp $SP --grid $GRID \
        --semicanon 0 --work pcc | tee step5_pcc.log
}

step6_triples() {
    local e_pcc
    e_pcc=$(python -c "import json; print(json.load(open('pcc/result.json'))['E_PCC'])")

    # coherent (T): projected (T) in the CC orbitals, effective-Fock orbital energies
    python $REPO/SUCCSD_new/run_pcc_t.py \
        --fc FC.p --mo pcc/MO_used.p --t1t2 pcc/T1T2.p \
        --newmo SGHFFock.p --moe eff \
        --sp $SP --grid $GRID --epcc $e_pcc --tag _coh | tee step6_coh_t.log

    # +(T): semicanonical orbitals of the bare Fock in the CC orbital spaces, then single-determinant (T)
    python $REPO/SUCCSD_new/semicanon_plusT.py FC.p pcc/MO_used.p $SP MO_sc.p | tee step6_semicanon.log
    python $REPO/SUCCSD_new/run_pcc_t.py \
        --fc FC.p --mo pcc/MO_used.p --t1t2 pcc/T1T2.p \
        --newmo MO_sc.p --moe fock \
        --sp 0 --grid 1 1 --epcc $e_pcc --tag _plus | tee step6_plus_t.log
}


# ----------------------------------------------------------------------------------------------- main
case "$1" in
    1) step1_frozen_core ;;
    2) step2_hf ;;
    3) step3_gradient_phf ;;
    4) step4_fock_phf ;;
    5) step5_pccsd ;;
    6) step6_triples ;;
    *) sed -n '2,15p' "$0"; exit 1 ;;
esac
