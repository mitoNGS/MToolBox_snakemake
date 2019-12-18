#!/usr/bin/env python
# -*- coding: UTF-8 -*-
from typing import List, Tuple


def get_mapping_regions(contigs_file) -> List[Tuple[int, int]]:
    """ TODO: fix docstring.

    Parameters
    ----------
    contigs_file
        Contig file with headings like:
            >Contig.1|1-10000
            >Contig.1|11000-15000

    Returns
    -------
    List of interval tuples:
        [(1,10000), (11000,15000)]
    """
    contigs = open(contigs_file, 'r')
    j = []
    l = contigs.readline()
    while l:
        if l.startswith('>'):
            y = l.strip().split('|')[1].split('-')
            j.append((int(y[0]), int(y[1])))
        l = contigs.readline()
    return j
