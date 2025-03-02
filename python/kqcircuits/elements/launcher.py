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


from kqcircuits.pya_resolver import pya
from kqcircuits.util.parameters import Param, pdt

from kqcircuits.elements.element import Element
from kqcircuits.defaults import default_layers


class Launcher(Element):
    """The PCell declaration for a launcher for connecting wirebonds.

    Default wirebond direction to west, waveguide to east. Uses default ratio a
    and b for scaling the gap.
    """

    s = Param(pdt.TypeDouble, "Pad width", 300, unit="μm")
    l = Param(pdt.TypeDouble, "Tapering length", 300, unit="μm")

    def produce_impl(self):
        # optical layer

        # keep the a/b ratio the same, but scale up a and b
        f = self.s / float(self.a)

        # shape for the inner conductor
        pts = [
            pya.DPoint(0, self.a / 2 + 0),
            pya.DPoint(self.l, f * (self.a / 2)),
            pya.DPoint(self.l + self.s, f * (self.a / 2)),
            pya.DPoint(self.l + self.s, -f * (self.a / 2)),
            pya.DPoint(self.l, -f * (self.a / 2)),
            pya.DPoint(0, -self.a / 2 + 0)
        ]

        shifts = [
            pya.DVector(0, self.b),
            pya.DVector(0, self.b * f),
            pya.DVector(self.b * f, self.b * f),
            pya.DVector(self.b * f, -self.b * f),
            pya.DVector(0, -self.b * f),
            pya.DVector(0, -self.b),
        ]
        pts2 = [p + s for p, s in zip(pts, shifts)]
        pts.reverse()
        shape = pya.DPolygon(pts + pts2)
        self.cell.shapes(self.get_layer("base_metal_gap_wo_grid")).insert(shape)

        # protection layer
        shifts = [
            pya.DVector(0, self.margin),
            pya.DVector(0, self.margin),
            pya.DVector(self.margin, self.margin),
            pya.DVector(self.margin, -self.margin),
            pya.DVector(0, -self.margin),
            pya.DVector(0, -self.margin),
        ]
        pts2 = [p + s for p, s in zip(pts2, shifts)]
        shape = pya.DPolygon(pts2)
        self.cell.shapes(self.get_layer("ground_grid_avoidance")).insert(shape)

        # add reference point
        self.add_port("", pya.DPoint(0, 0), pya.DVector(-1, 0))

        super().produce_impl()
