#!/bin/sh

style=$(printf %03d $1)

mkdir _tmp$style
cd    _tmp$style

n=1
for i in ../*-${style}.png
do
  counter=$(printf %06d $n)
  ln -s $i img${counter}.png
  n=$(($n+1))
done

ffmpeg -r 10 -sameq -i img%06d.png ../animation-${style}.mp4

rm -f img*png
cd ..
rmdir _tmp$style

