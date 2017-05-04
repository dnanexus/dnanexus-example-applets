#!/bin/bash -e

if [ -z "$1" ]
then
  echo Usage: update.sh applet_name
  exit
fi

p="${1%/}"

echo Updating "$p"...

if ! dx ls "Developer Applets": >/dev/null
then
  echo "Error running: dx ls \"Developer Applets\":"
  exit
fi

echo Clearing existing version...
( dx ls "Developer Applets":${p} >& /dev/null ) && dx rm "Developer Applets":${p}
( dx ls "Developer Applets":/sources/${p}.tar.gz >& /dev/null ) && dx rm "Developer Applets":/sources/${p}.tar.gz

echo Re-packaging...
rm -f "${p}.tar.gz"
tar zcvf "${p}.tar.gz" ${p}

echo Uploading...
dx upload "${p}.tar.gz" -p -o "Developer Applets":/sources/${p}.tar.gz

echo Re-building...
dx-build-applet ${p} -d "Developer Applets":${p}
