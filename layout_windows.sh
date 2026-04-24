#!/usr/bin/env bash

if ! $(command -v niri $>/dev/null); then
  echo "niri executable not found, bailing out"
  exit
fi

niri-set-window-width.sh 'Yorick 1' 630
sleep 0.1
niri-set-window-width.sh 'Yorick 2' 380
sleep 0.1
niri msg action focus-window --id $(niri-get-id-from-title.sh "Yorick 4")
niri msg action consume-or-expel-window-left
sleep 0.1
niri msg action focus-window --id $(niri-get-id-from-title.sh "DISPLAY=:0 GDK_SCALE ~/p/m/n/mavis-pray")
