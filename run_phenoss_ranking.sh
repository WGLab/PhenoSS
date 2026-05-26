#!/usr/bin/env bash

# ======================================================
# USER INPUTS (override via CLI if desired)
# ======================================================

INPUTFILE=""
OUTPUTFILE=""

MODE="oard_first"        # oard_only | oard_first | hpodb_first | hpodb_only
FREQ_ASSIGNMENT="extrinsic_ic"

METHOD="Resnik"
HP_DB_SQLITE="hp.db"
HPO_DB_PATH="./doc/database/hpo_frequency.csv"
GENE_CONVERSION=""
URL="https://rare.cohd.io/api"
DATASET_ID="2"

GENE_OF_INTEREST=""
GENE_OUTFILE=""

PY_SCRIPT="./doc/phenoss.py"   # python script name


# ======================================================
# ARGUMENT PARSER
# ======================================================

while [[ $# -gt 0 ]]; do
    case "$1" in
        --inputfile) INPUTFILE="$2"; shift 2 ;;
        --outputfile) OUTPUTFILE="$2"; shift 2 ;;
        --mode) MODE="$2"; shift 2 ;;
        --freq_assignment) FREQ_ASSIGNMENT="$2"; shift 2 ;;
        --method) METHOD="$2"; shift 2 ;;
        --hp_db_sqlite) HP_DB_SQLITE="$2"; shift 2 ;;
        --hpo_db_path) HPO_DB_PATH="$2"; shift 2 ;;
        --gene_conversion) GENE_CONVERSION=1; shift ;;   # <-- flag only
        --url) URL="$2"; shift 2 ;;
        --dataset_id) DATASET_ID="$2"; shift 2 ;;
        --gene_of_interest) GENE_OF_INTEREST="$2"; shift 2 ;;
        --gene_outfile) GENE_OUTFILE="$2"; shift 2 ;;
        *) echo "Unknown argument $1"; exit 1 ;;
    esac
done


# ======================================================
# CHECK REQUIRED
# ======================================================

if [[ -z "$INPUTFILE" || -z "$OUTPUTFILE" ]]; then
    echo "ERROR: --inputfile and --outputfile are required"
    exit 1
fi


# ======================================================
# RUN
# ======================================================

CMD="python ${PY_SCRIPT} \
    --inputfile ${INPUTFILE} \
    --outputfile ${OUTPUTFILE} \
    --mode ${MODE} \
    --freq_assignment ${FREQ_ASSIGNMENT} \
    --method ${METHOD} \
    --hp_db_sqlite ${HP_DB_SQLITE} \
    --hpo_db_path ${HPO_DB_PATH} \
    --url ${URL} \
    --dataset_id ${DATASET_ID}"

if [[ ! -z "$GENE_CONVERSION" ]]; then
    CMD="$CMD --gene_conversion"
fi

if [[ ! -z "$GENE_OF_INTEREST" ]]; then
    CMD="$CMD --gene_of_interest ${GENE_OF_INTEREST}"
fi

if [[ ! -z "$GENE_OUTFILE" ]]; then
    CMD="$CMD --gene_outfile ${GENE_OUTFILE}"
fi

echo "Running:"
echo $CMD
echo ""

eval $CMD
