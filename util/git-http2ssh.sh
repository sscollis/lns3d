#!/bin/bash
if [ $# -ne 2 ];
then
  echo Usage:  $0 username repo
  exit 1
fi
git remote -v
git remote set-url origin git@github.com:$1/$2.git
