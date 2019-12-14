#!/bin/bash

for i in {1..10}
do
  pdfcrop Presentation1-$i.pdf image$i.pdf 
done
