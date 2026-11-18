#!/bin/bash
# =====================================================================================================
#  N2 / 6-31G, R = 3.0 bohr, N 1s frozen
#  frozen core -> RHF / UHF -> SGHF -> gauge fixing -> SGCCSD -> coherent (T) and +(T)
# =====================================================================================================
#  Run one step at a time and check its log before going on:
#     bash run.sh 1    frozen core and integrals (Main_FC.py)               -> FC.p
#     bash run.sh 2    RHF / UHF in the frozen-core basis (Main_HF.py)      -> GHFMO.p (stable UHF)
#     bash run.sh 3    gradient SGHF, one basin hop, from the UHF           -> SGHF.p
#     bash run.sh 4    Fock-based SGHF from the gradient solution           -> SGHFFock.p (orbitals, eff.-Fock energies)
#     bash run.sh 5    gauge fixing (FixSz: <S_z> = 0), Fock SGHF again     -> SGHFFock.p in the fixed gauge
#     bash run.sh 6    SGCCSD in the effective-Fock orbitals                -> pcc/
#     bash run.sh 7    coherent (T) and +(T)                                -> pcc/result_T_coh.json, result_T_plus.json
#  Settings: parameter.py (molecule, frozen core, point group), parameter_fc.py (spin projection, grid 14 x 6).
#  The PCC module must be compiled first (README, "Compiling the code").
# =====================================================================================================
set -eo pipefail
cd "$(dirname "$0")"

REPO=$(cd ../.. && pwd)
export PYTHONPATH=$REPO
export LD_LIBRARY_PATH=$REPO/SUCCSD_new:$LD_LIBRARY_PATH
ln -sf $REPO/lebedev_*.txt .                  # SGHF: the Lebedev grids are read from the working folder

# spin projection of the PCC steps (as in parameter_fc.py)
SP=1
GRID="14 6"


# ----------------------------------------------------------------------------------------------- helper
# Run a script of the repository so that it imports parameter.py / parameter_fc.py of THIS folder
# ("python -c" puts the current folder first on the import path).
run_repo_script() {
    local script=$1
    shift
    python -c "import sys; sys.argv = ['$script'] + sys.argv[1:]; exec(open('$REPO/$script').read())" "$@"
}

# Fock-based PHF from the given guess; one retry from the lowest-energy iterate, then stop
fock_phf() {
    local guess=$1 log=$2
    run_repo_script PHFFock.py $guess | tee $log.log
    if grep -q "converged: True" $log.log; then
        return
    fi
    echo "Fock PHF not converged: retry once from the lowest-energy iterate"
    run_repo_script PHFFock.py newmo_best.p | tee ${log}_retry.log
    if ! grep -q "converged: True" ${log}_retry.log; then
        echo "STOP: no stationary PHF solution, do not run PCC on it"
        exit 1
    fi
}


# ----------------------------------------------------------------------------------------------- steps
step1_frozen_core() {
    run_repo_script Main_FC.py | tee step1_fc.log
}

step2_hf() {
    run_repo_script Main_HF.py | tee step2_hf.log
}

step3_gradient_phf() {
    # If this stops with "Not Hermitian": check that step 2 ended on a stable UHF, and use a larger
    # spin-projection grid in parameter_fc.py (e.g. ngrid = [26, 12]).
    rm -f SGHF.p hubsghf.p                    # restart files of an earlier PHF run
    python -c "
import PHF
PHF.nhop = 1           # one basin hop from the UHF (0: plain BFGS, stays at the nearest stationary point)
PHF.fixgauge = True    # SGHF: FixGauge (minimize the determinant <S^2> over spin gauges) at the end
exec(open('$REPO/Main_PHF.py').read())
" 2>&1 | tee step3_phf.log
    if [ ! -f SGHF.p ]; then
        # the Fortran code stops without an error status (e.g. "Not Hermitian" after a random basin-hop step)
        echo "STOP: the gradient PHF did not finish (see step3_phf.log); rerun this step"
        exit 1
    fi
}

step4_fock_phf() {
    fock_phf SGHF.p step4_fock
}

step5_gauge_fixing() {
    # FixSz (FixGauge.py): real boost exp(a S_z) with <S_z> = 0 of the determinant, then Fock PHF again (as for Cr2)
    python -c "
import pickle
import numpy as np
from parameter_fc import NAO, NOccSO
from Spin import calcS
from FixGauge import FixSz
data = pickle.load(open('SGHFFock.p', 'rb'))     # [MO, fsp, fpg, fk, moe]
mo = np.array(data[0], dtype=complex)
print('before FixSz: <Sx, Sy, Sz> =', np.real(calcS(mo, NOccSO, NAO)))
mo = FixSz(mo)
print('after FixSz:  <Sx, Sy, Sz> =', np.real(calcS(mo, NOccSO, NAO)))
data[0] = mo.real if abs(mo.imag).max() < 1e-10 else mo
pickle.dump(data, open('SGHF_sz.p', 'wb'))
" | tee step5_fixsz.log
    fock_phf SGHF_sz.p step5_fock
}

step6_pccsd() {
    # effective-Fock orbital energies (5th entry of SGHFFock.p) for the preconditioner
    python -c "import pickle, numpy; numpy.savetxt('moe.txt', pickle.load(open('SGHFFock.p', 'rb'))[4])"
    PCC_MOE_FILE=$PWD/moe.txt python $REPO/SUCCSD_new/run_pcc.py \
        --module PUCC_X2 --libdir $REPO/SUCCSD_new \
        --fc FC.p --mo SGHFFock.p --sp $SP --grid $GRID \
        --semicanon 0 --work pcc | tee step6_pcc.log
}

step7_triples() {
    local e_pcc
    e_pcc=$(python -c "import json; print(json.load(open('pcc/result.json'))['E_PCC'])")

    # coherent (T): projected (T) in the CC orbitals, effective-Fock orbital energies
    python $REPO/SUCCSD_new/run_pcc_t.py \
        --fc FC.p --mo pcc/MO_used.p --t1t2 pcc/T1T2.p \
        --newmo SGHFFock.p --moe eff \
        --sp $SP --grid $GRID --epcc $e_pcc --tag _coh | tee step7_coh_t.log

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
    5) step5_gauge_fixing ;;
    6) step6_pccsd ;;
    7) step7_triples ;;
    *) sed -n '2,16p' "$0"; exit 1 ;;
esac
