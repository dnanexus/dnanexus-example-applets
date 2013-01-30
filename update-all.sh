#!/bin/bash -e

for applet in *
do
  [ -d "$applet" ] && [ "$applet" != "Attic" ] && ./update.sh "$applet"
done
