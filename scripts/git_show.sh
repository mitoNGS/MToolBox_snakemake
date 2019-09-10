#!/usr/bin/env bash

cd $(dirname "$(readlink -f "$0")")
echo ""
echo "Hello! This is MToolBox. You are using this version:"
echo ""
echo $(git show --quiet --decorate | grep ^commit)
echo $(git show --quiet --decorate | grep ^Date)
echo ""
cd - > /dev/null