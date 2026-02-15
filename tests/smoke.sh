#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
BIN="$ROOT_DIR/dedup"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

fail() {
  echo "FAIL: $1" >&2
  exit 1
}

assert_eq() {
  local actual="$1"
  local expected="$2"
  local msg="$3"
  if [[ "$actual" != "$expected" ]]; then
    fail "$msg (expected=$expected actual=$actual)"
  fi
}

# 1) Similarity mode must catch near-duplicates even when first k bases differ.
cat > "$TMP_DIR/prefix_bug.fasta" << 'FASTA'
>s1
ACGTACGTACGTACGTACGT
>s2
TCGTACGTACGTACGTACGT
FASTA
"$BIN" -i "$TMP_DIR/prefix_bug.fasta" -o "$TMP_DIR/prefix_bug_out.fa" -p 95 -t 2 >/dev/null 2>&1
assert_eq "$(grep -c '^>' "$TMP_DIR/prefix_bug_out.fa")" "1" "prefix mismatch duplicate should be removed"

# 2) Empty sequence duplicates should be removed in similarity mode.
cat > "$TMP_DIR/empty.fasta" << 'FASTA'
>e1

>e2

FASTA
"$BIN" -i "$TMP_DIR/empty.fasta" -o "$TMP_DIR/empty_out.fa" -p 95 -t 2 >/dev/null 2>&1
assert_eq "$(grep -c '^>' "$TMP_DIR/empty_out.fa")" "1" "empty duplicate should be removed"

# 3) Keep first occurrence deterministically in exact mode.
cat > "$TMP_DIR/order.fasta" << 'FASTA'
>first
ACGTACGTACGTACGT
>second
ACGTACGTACGTACGT
FASTA
for _ in $(seq 1 10); do
  "$BIN" -i "$TMP_DIR/order.fasta" -o "$TMP_DIR/order_out.fa" -p 100 -t 8 >/dev/null 2>&1
  assert_eq "$(head -n 1 "$TMP_DIR/order_out.fa")" ">first" "exact mode should keep first occurrence"
done

# 4) Unique inputs should all be retained.
{
  for i in $(seq 0 199); do
    echo ">u$i"
    printf 'ACGT%06dTGCA\n' "$i"
  done
} > "$TMP_DIR/unique.fasta"
"$BIN" -i "$TMP_DIR/unique.fasta" -o "$TMP_DIR/unique_out.fa" -p 100 -t 8 >/dev/null 2>&1
assert_eq "$(grep -c '^>' "$TMP_DIR/unique_out.fa")" "200" "unique inputs should all be kept"

echo "smoke tests passed"
