#!/bin/bash
# Starts the quicklooks HTTP server on port 8765 if it is not already running.
# Intended to be called from cron, e.g.:
#   */5 * * * * /home/cwalden/git/kepler-radar-utils/start_quicklooks_server.sh

QUICKLOOKS_DIR="/data/processing/kepler/reading-general/quicklooks"
PORT=8765

pgrep -f "http.server ${PORT}" > /dev/null && exit 0

cd "${QUICKLOOKS_DIR}" || exit 1
nohup python -m http.server ${PORT} > /dev/null 2>&1 &
