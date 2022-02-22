#!/bin/bash
for i in {5..175}
do
  ./ExportForTagging ${i} ${i} > m12b/${i}.csv
done
