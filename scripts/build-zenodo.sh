#!/usr/bin/env bash
# Build the Zenodo deposition bundle in _zenodo/.
#
#   bash scripts/build-zenodo.sh
#
# SCOPE. West Coast Arborists agreed to release a subset of the inventory with
# per-tree location removed. Only the public, location-free release is deposited.
# The private preprocessing intermediates (wca_filtered.csv,
# wca_filtered_annotated.csv, wca-partitions.RData) all carry latitude/longitude
# and MUST NOT be included -- an earlier version of this script bundled them, and
# the check at the end exists so that cannot happen again silently.
#
# Everything here is derived from files already in the working tree; the script is
# idempotent and safe to re-run.
set -euo pipefail
cd "$(dirname "$0")/.."

OUT="_zenodo"
USFS_ROOT="data/_RDS-2016-0005"
USFS="$USFS_ROOT/Data"

# built by scripts/py/3_make_public_data.py
PUBLIC_FILES=(wca_inventory_public.csv species_frost_free_days.csv DATA_DICTIONARY.md)
# the three the analysis reads; the whole archive ships regardless
USFS_FILES=(TS1_Regional_information.csv TS3_Raw_tree_data.csv TS6_Growth_coefficients.csv)

mkdir -p "$OUT"

# the dictionary is authored documentation and is tracked at the repo root
# (docs/ is gitignored for pkgdown); copy it in at build time so there is one copy
cp PUBLIC_DATA.md "$OUT/DATA_DICTIONARY.md"

for f in "${PUBLIC_FILES[@]}"; do
  [ -f "$OUT/$f" ] || { echo "missing: $OUT/$f -- run scripts/py/3_make_public_data.py" >&2; exit 1; }
done
for f in "${USFS_FILES[@]}"; do
  [ -f "$USFS/$f" ] || { echo "missing: $USFS/$f" >&2; exit 1; }
done

echo "building $OUT/usfs-rds-2016-0005.zip ..."
# ship RDS-2016-0005 whole and unmodified: all nine TS tables plus the official
# metadata and the Jan 2020 errata, so the deposition carries the citable unit.
usfs_tmp="$(mktemp -d)"
cp -Rp "$USFS_ROOT" "$usfs_tmp/RDS-2016-0005"
find "$usfs_tmp" -name '.DS_Store' -delete
( cd "$usfs_tmp" && zip -q -r -X out.zip RDS-2016-0005 )
mv "$usfs_tmp/out.zip" "$OUT/usfs-rds-2016-0005.zip"
rm -rf "$usfs_tmp"

sum()  { md5 -q "$1" 2>/dev/null || md5sum "$1" | cut -d' ' -f1; }
hsize(){
  local b; b=$(stat -f %z "$1" 2>/dev/null || stat -c %s "$1")
  awk -v b="$b" 'BEGIN{ split("B KB MB GB TB", u, " ")
    i=1; while (b >= 1024 && i < 5) { b/=1024; i++ }
    printf (i==1 ? "%d %s" : "%.1f %s"), b, u[i] }'
}

{
cat <<'HEADER'
================================================================================
Supporting data for "New allometric equations for urban trees"
https://github.com/ruizt/allometry
================================================================================

FILES

HEADER

printf '  %-44s  %11s  %s\n' "file" "size" "md5"
printf '  %-44s  %11s  %s\n' "--------------------------------------------" "-----------" "--------------------------------"
for f in "${PUBLIC_FILES[@]}"; do
  printf '  %-44s  %11s  %s\n' "$f" "$(hsize "$OUT/$f")" "$(sum "$OUT/$f")"
done
printf '  %-44s  %11s  %s\n' "usfs-rds-2016-0005.zip" "$(hsize "$OUT/usfs-rds-2016-0005.zip")" "$(sum "$OUT/usfs-rds-2016-0005.zip")"

cat <<'FILENOTES'

  wca_inventory_public.csv     Tree measurements from the California urban tree
                               inventory: species, DBH, height, date, estimated
                               age, canopy width, and per-species outlier flags.
                               One row per measurement; repeated measurements of a
                               tree share a tree_id. No coordinates.

  species_frost_free_days.csv  Mean annual frost-free days per species.

  DATA_DICTIONARY.md           Column definitions for both CSVs.

  usfs-rds-2016-0005.zip       The USDA Forest Service Urban Tree Database and
                               Allometric Equations archive, complete and
                               unmodified: nine data tables, official metadata,
                               and the January 2020 errata.

PROVENANCE

  Tree records were contributed by West Coast Arborists, Inc. (WCA) from their
  California urban tree inventory, and are released here by agreement.

  Canopy widths are remotely sensed from USDA National Agriculture Imagery
  Program (NAIP) imagery, 2020 and 2022, and are not WCA data.

  Frost-free days are derived from the USDA soil frost-free days raster.

  Ages and outlier flags were computed by the preprocessing code in the
  repository above.

  The USFS archive is redistributed unchanged for convenience. Cite the original:
    McPherson, van Doorn & Peper. Urban Tree Database and Allometric Equations.
    Forest Service Research Data Archive, RDS-2016-0005.
    https://doi.org/10.2737/RDS-2016-0005

FILENOTES

echo "Generated $(date '+%Y-%m-%d') by scripts/build-zenodo.sh"
} > "$OUT/description.txt"

# ---- safety check: nothing in the deposition may carry coordinates ----------
echo
echo "checking the deposition for location data ..."
leak=0
for f in "$OUT"/*.csv; do
  [ -e "$f" ] || continue
  if head -1 "$f" | tr ',' '\n' | grep -qiE '^"?(lat|lon)'; then
    echo "  LEAK: $f has a location column" >&2; leak=1
  fi
done
for f in "$OUT"/*.zip; do
  [ -e "$f" ] || continue
  if unzip -l "$f" | grep -qiE 'wca_filtered|wca-partitions|wca_name_resolved|wca_removals'; then
    echo "  LEAK: $f contains a private inventory product" >&2; leak=1
  fi
done
[ "$leak" -eq 0 ] && echo "  ok: no location columns, no private inventory products"
[ "$leak" -eq 0 ] || { echo "REFUSING to finish with a leak staged" >&2; exit 1; }

echo
echo "done:"
ls -la "$OUT"
