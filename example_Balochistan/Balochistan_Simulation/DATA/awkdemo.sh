#! /bin/bash
awk '{if($1~/S/){print $1 $2 $3}}' Par_file_faults
