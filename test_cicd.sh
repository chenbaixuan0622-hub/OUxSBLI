#!/bin/bash
# Compile-test all CICD build targets.
# Run from the repository root: bash test_cicd.sh
# Requires: NVIDIA HPC SDK with mpif90 on PATH.

PASS=()
FAIL=()

run_case() {
    local dir="$1"
    echo "=== Building $dir ==="
    if (cd "$dir" && make clean -s && make 2>&1); then
        PASS+=("$dir")
        echo "--- PASS: $dir ---"
    else
        FAIL+=("$dir")
        echo "--- FAIL: $dir ---"
    fi
    echo
}

# Cartesian solver: 9 cases (KEEP/SLAU/Hybrid × 2nd/4th/6th order)
for case in CICD_KEEP2 CICD_KEEP4 CICD_KEEP6 \
            CICD_SLAU2 CICD_SLAU4 CICD_SLAU6 \
            CICD_Hybrid2 CICD_Hybrid4 CICD_Hybrid6; do
    run_case "3D_solver/CICD/$case"
done

# Curvilinear solver: 3 cases (KEEP/SLAU/Hybrid × 2nd order only)
for case in CICD_KEEP2 CICD_SLAU2 CICD_Hybrid2; do
    run_case "3D_solver_curv/CICD/$case"
done

echo "=============================="
echo "PASSED (${#PASS[@]}):"
for d in "${PASS[@]}"; do echo "  $d"; done
echo "FAILED (${#FAIL[@]}):"
for d in "${FAIL[@]}"; do echo "  $d"; done
echo "=============================="

[ ${#FAIL[@]} -eq 0 ]
