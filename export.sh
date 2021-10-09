#!/bin/bash
for i in {9..255}
do
  ./ExportForTagging ${i} ${i} > ${i}.csv
done
