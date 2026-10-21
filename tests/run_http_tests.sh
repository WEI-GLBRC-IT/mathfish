#!/bin/bash
# POST every analysis form through Apache + PHP, using the same inputs as
# run_cli_tests.sh (which calls the Python scripts directly and so never
# exercises PHP). Intended to be invoked with:
#   docker cp tests/run_http_tests.sh mathfish:/tmp/run_http_tests.sh
#   docker exec mathfish bash /tmp/run_http_tests.sh
#
# A page passes when it returns HTTP 200, shows no error text, and every plot
# it links to downloads as a PNG. PHP warnings go to the Apache error log, not
# the page, so also check `docker logs` (CI does).

set -o pipefail

BASE="${BASE:-http://localhost}"

ORG='AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGCAGCTTGCTGCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTAGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACACG'
# 2-mismatch variant (positions 346-347 AC -> AA in probe binding region)
MM='AAATTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAGCAGCTTGCTGCTTTGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTAGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTAGTAGGTGGGGTAACGGCTCACCTAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTAAGGCAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACACG'

PROBE='GCTGCCTCCCGTAGGAGTCTGGAC'
COMPE='GCAGCCACACTGGAACTGAGACAC'

checked=0
passed=0

# check_page <name> <page.php> <curl --data-urlencode args...>
check_page() {
    local name="$1"; shift
    local page="$1"; shift
    local body code fail="" plots=0 src got
    checked=$((checked + 1))
    body=$(mktemp)

    code=$(curl -sS -o "$body" -w '%{http_code}' "$@" "$BASE/$page") || fail="curl failed"
    [ "$code" = 200 ] || fail="${fail:+$fail; }HTTP $code"
    if grep -qE 'Error encountered|Long sequence error|Image not found|<b>(Warning|Fatal error|Deprecated|Notice)</b>' "$body"; then
        fail="${fail:+$fail; }error text in page: $(grep -oE 'Error encountered|Long sequence error|Image not found|<b>(Warning|Fatal error|Deprecated|Notice)</b>' "$body" | head -1)"
    fi
    for src in $(grep -oE 'src="userFiles/[A-Za-z0-9_]+\.png"' "$body" | cut -d'"' -f2); do
        plots=$((plots + 1))
        got=$(curl -sS -o /dev/null -w '%{http_code} %{content_type}' "$BASE/$src")
        [ "$got" = "200 image/png" ] || fail="${fail:+$fail; }$src -> $got"
    done
    [ "$plots" -gt 0 ] || fail="${fail:+$fail; }no plot linked"
    rm -f "$body"

    if [ -z "$fail" ]; then
        passed=$((passed + 1))
        echo "[$name] OK ($page, $plots plot(s))"
    else
        echo "[$name] FAIL ($page): $fail"
    fi
}

common=(--data-urlencode T=46 --data-urlencode NA=1)
target=(--data-urlencode targetmolecule=ssu --data-urlencode domain=b)
series=(--data-urlencode "organism=$ORG" --data-urlencode ts1=196 --data-urlencode ts2=213 --data-urlencode pl=15)

check_page general Step2.php "${common[@]}" "${target[@]}" --data-urlencode Po=250 --data-urlencode entry=p \
    --data-urlencode "ipprobesequence=$PROBE" --data-urlencode "organism=$ORG"
check_page dg1 dg1.php "${common[@]}" "${series[@]}"
check_page dg2 dg2.php "${common[@]}" "${series[@]}"
check_page dg3 dg3.php "${common[@]}" "${series[@]}" "${target[@]}"
check_page dgov dgovser.php "${common[@]}" "${series[@]}" "${target[@]}"
check_page formamide formamide.php "${common[@]}" "${target[@]}" --data-urlencode Po=250 \
    --data-urlencode "probe[]=$PROBE" --data-urlencode "organism[]=$ORG" --data-urlencode chkdata=0
check_page mismatch Step3.php "${common[@]}" "${target[@]}" --data-urlencode Po=250 \
    --data-urlencode "ipprobesequence=$PROBE" --data-urlencode "organism=$ORG" --data-urlencode "mismatch=$MM"
check_page competitor compe.php "${common[@]}" "${target[@]}" --data-urlencode Po=250 --data-urlencode Co=5000 \
    --data-urlencode "probe=$PROBE" --data-urlencode "compe=$COMPE" \
    --data-urlencode "organism=$ORG" --data-urlencode "nontarget=$ORG"

echo "pages checked: $checked, passed: $passed"
[ "$checked" -eq 8 ] && [ "$passed" -eq 8 ]
