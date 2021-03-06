#!/usr/bin/env python
# -*- coding: UTF-8 -*-
from modules.classifier import tree


# TODO: rm_file is not used anywhere
def write_haplogroup(file_handle, rm_file, haplogroup):
    # solo quelli filtrati
    filtered_positions = tree.filter_positions(haplogroup)
    print("CIAO")
    print(haplogroup)
    print(filtered_positions)
    pos_list = []
    for position in filtered_positions:
        pos_list.append((haplogroup.name,) + position.print_table())
    file_handle.writerows(pos_list)
