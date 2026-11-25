#!/usr/bin/env bash
# Cell-Hub end-to-end test of a Docker image: database checks + every module driven in a headless browser on the demo data.
# Usage: bash run_cellhub_tests.sh [image] [demo_data_dir] [modules]
#   image          default gaspardmacaux/cell-hub:2.0.0
#   demo_data_dir  default ~/Bureau/cellhub_demo/cellhub_demo_data (files of the Zenodo record)
#   modules        default single,cellchat,trajectory,multiple,spatial
# Output: cellhub_test_<date>/report.html (+ screenshots/, downloads/, logs)
set -uo pipefail
IMAGE="${1:-gaspardmacaux/cell-hub:2.0.0}"
DEMO="$(realpath "${2:-$HOME/Bureau/cellhub_demo/cellhub_demo_data}")"
MODULES="${3:-single,cellchat,trajectory,multiple,spatial}"
# APP_DIR=/path/to/mon_app tests that app folder (files with the _docker names, databases/ included) inside the image,
# mounted over the app of the image: the code can be checked before a rebuild, with the packages of the image
APP_DIR="${APP_DIR:-}"
APP_MOUNT=()
if [ -n "$APP_DIR" ]; then
  APP_DIR="$(realpath "$APP_DIR")"
  [ -d "$APP_DIR/databases" ] || echo "WARNING: $APP_DIR has no databases/ folder: CellChat and the database checks will fail"
  APP_MOUNT=(-v "$APP_DIR":/srv/shiny-server/mon_app:ro)
fi
TESTS="$(cd "$(dirname "$0")" && pwd)"
OUT="$PWD/cellhub_test_$(date +%Y%m%d_%H%M)"
PW_IMAGE="mcr.microsoft.com/playwright/python:v1.48.0-jammy"
APP="cellhub_test_app"; NET="cellhub_test_net"
mkdir -p "$OUT"
for f in pbmc_1k_v3_filtered_feature_bc_matrix.h5 pbmc_1k_v2_filtered_feature_bc_matrix.h5 visium_mouse_brain.zip; do
  [ -f "$DEMO/$f" ] || { echo "ERROR: $DEMO/$f not found (demo data folder)"; exit 1; }
done
echo "Image: $IMAGE | demo data: $DEMO | modules: $MODULES | output: $OUT"
docker image inspect "$IMAGE" --format 'Image id: {{.Id}} created {{.Created}}' > "$OUT/image.txt" 2>/dev/null || { echo "ERROR: image $IMAGE not found locally (build it or check the tag with: docker images)"; exit 1; }
cat "$OUT/image.txt"
# Browser image + Python playwright package, built once on the host network (containers on custom networks may lack DNS)
if ! docker image inspect cellhub-uitest:1.48.0 >/dev/null 2>&1; then
  printf 'FROM mcr.microsoft.com/playwright/python:v1.48.0-jammy\nRUN pip install --no-cache-dir playwright==1.48.0\n' | docker build --network host -q -t cellhub-uitest:1.48.0 - || { echo "ERROR: could not build the test browser image"; exit 1; }
fi
PW_IMAGE="cellhub-uitest:1.48.0"

# 1. Database checks with the app's own helpers (independent of the interface)
echo "== Database checks"
docker run --rm --entrypoint Rscript --user "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "$TESTS":/tests:ro -v "$OUT":/out "${APP_MOUNT[@]}" "$IMAGE" /tests/db_test.R /srv/shiny-server/mon_app /out 2>&1 | tee "$OUT/db_test.log" | tail -3

# 1b. Functional test: the app's functions on the demo data, compared with plain Seurat / CellChat / harmony (FUNCTIONAL=0 skips it)
if [ "${FUNCTIONAL:-1}" = "1" ]; then
  echo "== Functional test (app functions vs independent reference, ~20-40 min)"
  docker run --rm --entrypoint Rscript --user "$(id -u):$(id -g)" -e HOME=/tmp \
    -v "$TESTS":/tests:ro -v "$DEMO":/demo:ro -v "$OUT":/out "${APP_MOUNT[@]}" "$IMAGE" /tests/functional_test.R /srv/shiny-server/mon_app /demo /out 2>&1 \
    | tee "$OUT/functional_test.log" | grep -E "^\[(OK|KO|WARN|SKIP)\]|^Functional test|====" 
fi

# 2. App container on a private network (also visible at http://localhost:3839 to watch it)
echo "== Starting the app"
docker rm -f "$APP" >/dev/null 2>&1
docker network create "$NET" >/dev/null 2>&1
docker run -d --name "$APP" --network "$NET" -p 3839:3838 "${APP_MOUNT[@]}" "$IMAGE" >/dev/null
for i in $(seq 1 120); do curl -sf -o /dev/null http://localhost:3839 && break; sleep 5; done
curl -sf -o /dev/null http://localhost:3839 || { echo "ERROR: the app did not start"; docker logs "$APP" > "$OUT/app_log.txt" 2>&1; exit 1; }

# 3. Memory of the app container every 5 s (epoch, GB), matched to the steps in the report
( while docker inspect "$APP" >/dev/null 2>&1; do
    docker stats --no-stream --format '{{.MemUsage}}' "$APP" 2>/dev/null | awk -v t="$(date +%s)" '{
      v=$1; u=v; gsub(/[0-9.]/,"",u); gsub(/[A-Za-z]/,"",v);
      f=(u=="GiB")?1:(u=="MiB")?1/1024:(u=="KiB")?1/1048576:(u=="TiB")?1024:0; printf "%s,%.3f\n", t, v*f }' >> "$OUT/mem.csv"
    sleep 5
  done ) &
MEM_PID=$!

# 4. Interface: every module in a headless browser
echo "== Interface tests (this takes a while; progress below)"
# Runs as root inside the test browser image; ownership of the results restored after
docker run --rm --network "$NET" --ipc=host \
  -v "$TESTS":/tests:ro -v "$DEMO":/demo:ro -v "$OUT":/out "$PW_IMAGE" \
  bash -c "python /tests/ui_test.py --url http://$APP:3838 --demo /demo --out /out --modules $MODULES; chown -R $(id -u):$(id -g) /out" 2>&1 | tee "$OUT/ui_test.log"

# 5. Logs, cleanup and report
docker logs "$APP" > "$OUT/app_log.txt" 2>&1
docker rm -f "$APP" >/dev/null 2>&1; kill "$MEM_PID" 2>/dev/null; docker network rm "$NET" >/dev/null 2>&1
docker run --rm --user "$(id -u):$(id -g)" -v "$TESTS":/tests:ro -v "$OUT":/out "$PW_IMAGE" python /tests/report.py /out
echo "Open: $OUT/report.html"
