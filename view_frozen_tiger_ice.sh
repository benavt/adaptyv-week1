#!/usr/bin/env bash
cd "$(dirname "$0")"
pymol "01_Staging/frozen-tiger-ice/frozen_tiger_ice_top5.pse" \
  -d "load 02_Design/foundry/runs/DD2_frozen-tiger-ice_hotspots_09_29/run_display/frozen_tiger_ice_hotspots.pse, partial=1"
