#!/bin/bash

BRANCH=`git branch 2> /dev/null | grep "*" | sed 's#*\ \(.*\)#\1#'`
if [ "$BRANCH" != "master" ]
then
    exit 1
fi

DESC=`git describe`

TMP=`mktemp -d temp.XXXXXX.$$`

cp -r ./src/biograpy/docs/_build/* $TMP 
#rm -f $TMP/*~
git checkout gh-pages
cp -r $TMP/* ./
git add -A
git commit -m "Update docs to $DESC"
git push origin gh-pages
git checkout master
rm -rf $TMP
