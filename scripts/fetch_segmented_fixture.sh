#!/usr/bin/env bash
# Freeze a GenBank XML snapshot for the segmented (influenza A, taxid 11320)
# test, so `segmented_xml_test` never touches NCBI.
#
# Why this exists: on a fresh run FETCH_GENBANK is 37.8% of segmented_test's
# total task time (129s of 341s measured), and it is the only step whose output
# depends on what NCBI returns on the day. Both matter - the second more than
# the first, because a change in NCBI's response format has already produced a
# CI failure that could not be reproduced locally.
#
# Refresh deliberately, not on a schedule: the point of a fixture is that it
# does not move under you.
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TAX_ID="11320"
OUT_DIR="${PROJECT_DIR}/test_data/iav_11320"
XML_SUBDIR="GenBank-XML"
REF_LIST="${PROJECT_DIR}/generic/influenza/ref_list_refmast.txt"
BATCH_SIZE="200"
EMAIL="${VGTK_EMAIL:-}"

usage() {
    cat <<EOF
Usage: $(basename "$0") [options]   (influenza A, taxid 11320)

Options:
  --email EMAIL        Email for NCBI E-utilities (or set VGTK_EMAIL)
  --out-dir PATH       Fixture root directory (default: test_data/iav_11320)
  --batch-size N       GenBankFetcher batch size (default: 200)
  --clean              Remove existing XML snapshot before downloading
  -h, --help           Show this help
EOF
}

CLEAN=0
while [[ $# -gt 0 ]]; do
    case "$1" in
        --email)
            EMAIL="$2"
            shift 2
            ;;
        --out-dir)
            OUT_DIR="$2"
            shift 2
            ;;
        --batch-size)
            BATCH_SIZE="$2"
            shift 2
            ;;
        --clean)
            CLEAN=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown argument: $1" >&2
            usage
            exit 1
            ;;
    esac
done

if [[ -z "${EMAIL}" ]]; then
    echo "--email (or VGTK_EMAIL) is required" >&2
    exit 1
fi

XML_DIR="${OUT_DIR}/${XML_SUBDIR}"
if [[ "${CLEAN}" -eq 1 ]]; then
    rm -rf "${XML_DIR}"
fi

mkdir -p "${OUT_DIR}"

python "${PROJECT_DIR}/scripts/GenBankFetcher.py" \
    --taxid "${TAX_ID}" \
    -o "${OUT_DIR}" \
    -d "${XML_SUBDIR}" \
    -b "${BATCH_SIZE}" \
    -e "${EMAIL}" \
    --ref_list "${REF_LIST}"

if [[ ! -d "${XML_DIR}" ]]; then
    echo "Expected XML directory not found: ${XML_DIR}" >&2
    exit 1
fi

XML_COUNT="$(find "${XML_DIR}" -maxdepth 1 -type f -name '*.xml' | wc -l | tr -d ' ')"
MANIFEST="${OUT_DIR}/manifest.tsv"
CHECKSUMS="${OUT_DIR}/checksums.sha256"

{
    echo -e "key\tvalue"
    echo -e "tax_id\t${TAX_ID}"
    echo -e "organism\tH10N8 subtype"
    echo -e "downloaded_utc\t$(date -u +'%Y-%m-%dT%H:%M:%SZ')"
    echo -e "xml_directory\t${XML_SUBDIR}"
    echo -e "xml_file_count\t${XML_COUNT}"
    echo -e "ref_list\t${REF_LIST}"
} > "${MANIFEST}"

if [[ "${XML_COUNT}" -gt 0 ]]; then
    find "${XML_DIR}" -maxdepth 1 -type f -name '*.xml' -print0 \
        | sort -z \
        | xargs -0 sha256sum > "${CHECKSUMS}"
else
    : > "${CHECKSUMS}"
fi

echo "Fixture ready: ${OUT_DIR}"
echo "XML files: ${XML_COUNT}"