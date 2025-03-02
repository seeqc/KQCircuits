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
from autologging import logged, traced

from kqcircuits.elements.airbridges.airbridge import Airbridge


@traced
@logged
class AirbridgeRectangular(Airbridge):
    """PCell declaration for a rectangular airbridge.

    Origin is at the geometric center. The airbridge is in vertical direction.

    Bottom parts of pads in bottom layer, bridge and top parts of pads in top layer. Pads and bridge are rectangular.
    Refpoints "port_a" and "port_b" at top pad points closest to origin.
    """

    default_type = "Airbridge Rectangular"

    bridge_width = Param(pdt.TypeDouble, "Bridge width", 20, unit="μm")

    def produce_impl(self):
        # shorthand
        (w, h, l, b, e) = (self.pad_width, self.pad_length, self.bridge_length, self.bridge_width, self.pad_extra)

        pts = [
            pya.DPoint(-self.pad_width / 2 - e, l / 2 - e),
            pya.DPoint(-self.pad_width / 2 - e, h + l / 2 + e),
            pya.DPoint(self.pad_width / 2 + e, h + l / 2 + e),
            pya.DPoint(self.pad_width / 2 + e, l / 2 - e),
        ]
        self._produce_bottom_pads(pts)

        # top layer
        pts = [
            pya.DPoint(-w / 2, h + l / 2),
            pya.DPoint(w / 2, h + l / 2),
            pya.DPoint(w / 2, l / 2),
            pya.DPoint(b / 2, l / 2),
            pya.DPoint(b / 2, -l / 2),
            pya.DPoint(w / 2, -l / 2),
            pya.DPoint(w / 2, -h - l / 2),
            pya.DPoint(-w / 2, -h - l / 2),
            pya.DPoint(-w / 2, -l / 2),
            pya.DPoint(-b / 2, -l / 2),
            pya.DPoint(-b / 2, l / 2),
            pya.DPoint(-w / 2, l / 2),
        ]
        self._produce_top_pads_and_bridge(pts)
