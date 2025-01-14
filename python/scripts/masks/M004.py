# This code is part of KQCircuits
# Copyright (C) 2021 IQM Finland Oy
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see
# https://www.gnu.org/licenses/gpl-3.0.html.
#
# The software distribution should follow IQM trademark policy for open-source software
# (meetiqm.com/developers/osstmpolicy). IQM welcomes contributions to the code. Please see our contribution agreements
# for individuals (meetiqm.com/developers/clas/individual) and organizations (meetiqm.com/developers/clas/organization).


from kqcircuits.chips.airbridge_crossings import AirbridgeCrossings
from kqcircuits.chips.quality_factor import QualityFactor
from kqcircuits.defaults import TMP_PATH
from kqcircuits.klayout_view import KLayoutView
from kqcircuits.masks.mask_set import MaskSet

"""M004 mask.

Q and AB tests

"""

view = KLayoutView(current=True, initialize=True)
layout = KLayoutView.get_active_layout()

m004 = MaskSet(layout, name="M004", version=2, with_grid=False)

box_map = {"A": [
    ["AB1", "AB2", "QSG"],
    ["QSA", "QSC", "QDG"],
    ["QDA", "QDC", "QDD"],
]}

mask_map = [
    ["A", "A", "A", "A", "A"],
    ["A", "A", "A", "A", "A"],
    ["A", "A", "A", "A", "A"],
    ["A", "A", "A", "A", "A"],
    ["A", "A", "A", "A", "A"],
]

m004.add_mask_layout(MaskSet.chips_map_from_box_map(box_map, mask_map))

parameters_qd = {
    "res_lengths": [4649.6, 4743.3, 4869.9, 4962.9, 5050.7, 5138.7, 5139., 5257., 5397.4, 5516.8, 5626.6, 5736.2,
                    5742.9, 5888.7, 6058.3, 6202.5, 6350., 6489.4],
    "type_coupler": ["interdigital", "interdigital", "interdigital", "gap", "gap", "gap", "interdigital", "interdigital", "interdigital", "gap",
                     "gap", "gap", "interdigital", "interdigital", "interdigital", "interdigital", "gap", "gap"],
    "l_fingers": [19.9, 54.6, 6.7, 9.7, 22.8, 30.5, 26.1, 14.2, 18.2, 10.9, 19.8, 26.4, 34.2, 19.9, 25.3, 8., 15.8,
                  22.2],
    "n_fingers": [4, 2, 2, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 2, 2, 4, 4],
    "res_beg": ["galvanic"]*18
}

parameters_qs = {
    "res_lengths": [4649.6, 4908.9, 5208.5, 5516.8, 5848.9, 6217.4],
    "type_coupler": ["interdigital", "interdigital", "interdigital", "gap", "gap", "gap"],
    "l_fingers": [19.9, 7.3, 15.2, 10.9, 18.5, 23.6],
    "n_fingers": [4, 4, 2, 4, 4, 4],
    "res_beg": ["galvanic"]*6
}

m004.add_chip(AirbridgeCrossings, "AB1", crossings=1)
m004.add_chip(AirbridgeCrossings, "AB2", crossings=10)
m004.add_chip(QualityFactor, "QSG", **parameters_qs, n_ab=6*[0], res_term=6*["galvanic"])
m004.add_chip(QualityFactor, "QSA", **parameters_qs, n_ab=6*[0], res_term=6*["airbridge"])
m004.add_chip(QualityFactor, "QSC", **parameters_qs, n_ab=6*[5], res_term=6*["galvanic"])
m004.add_chip(QualityFactor, "QDG", **parameters_qd, n_ab=18*[0], res_term=18*["galvanic"])
m004.add_chip(QualityFactor, "QDA", **parameters_qd, n_ab=18*[0], res_term=18*["airbridge"])
m004.add_chip(QualityFactor, "QDC", **parameters_qd, n_ab=18*[5], res_term=18*["galvanic"])
m004.add_chip(QualityFactor, "QDD", **parameters_qd, n_ab=18*[5], res_term=18*["galvanic"])

m004.build()
m004.export(TMP_PATH, view)