#!/usr/bin/env bash
# Created by: Li Yao (ly349@cornell.edu)
# Created on: 2020-05-27
#
# Part of the PINTS project
# Copyright (C) 2020 Li Yao at the Yu Lab
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# This wrapper restricts Tfit to create only 16 threads instead of using all cores
set -o errexit;
export OMP_NUM_THREADS=16;
/local/storage/ly349/toolshed/BioQueue/workspace/5/bin/Tfit bidir -i $1 -j $2 -o . -N $3;
/local/storage/ly349/toolshed/BioQueue/workspace/5/bin/Tfit model -i $1 -j $2 -k *prelim_bidir_hits.bed -o . -N $3;