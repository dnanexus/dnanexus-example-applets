#!/bin/bash -e

for applet in *
do
  [ -d "$applet" ] && [ "$applet" != "Attic" ] && [ "$applet" != "bedtools" ] && ./update.sh "$applet"
done
