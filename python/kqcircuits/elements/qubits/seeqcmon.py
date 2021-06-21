import numpy as np
from math import pi as π
from autologging import traced

from kqcircuits.util.parameters import Param, pdt
from kqcircuits.elements.qubits.qubit import Qubit
from kqcircuits.pya_resolver import pya


@traced
class Seeqcmon(Qubit):
    """The PCell declaration for a Seeqc qubit.
    """

    body_radius = Param(
        pdt.TypeDouble, "Qubit body radius (μm)", 60, unit='μm')
    body_ring_gap = Param(pdt.TypeDouble, "Body-Ring gap (μm)", 5, unit='μm')
    ring_gnd_gap = Param(pdt.TypeDouble, "Ring-Gnd gap (μm)", 30, unit='μm')
    num_coupling = Param(pdt.TypeInt, "Number of couplings", 3)
    coupling_ring_width = Param(
        pdt.TypeDouble, "Coupling ring width (μm)", 5, unit='μm')
    stub_width = Param(pdt.TypeDouble, "Stub width (μm)", 2, unit='μm')
    stub_angle = Param(pdt.TypeDouble, "Stub angle (deg)", 10)

    def produce_impl(self):

        self._produce_body()

        self._produce_ring()

        return super().produce_impl()

    def _circle_reg_from_radius(self, r: float = 10, n_pts: int = 128):

        pts = [pya.DPoint(r*np.cos(φ), r*np.sin(φ))
               for φ in np.linspace(0, 2*π, n_pts)]
        poly = pya.DPolygon(pts)
        reg = pya.Region([poly.to_itype(self.layout.dbu)])

        return reg

    def _produce_arc_pts(self, r: float = 10, n_pts: int = 128, φ: float = 90, φ0=0) -> list:
        """Produce points for arc

        :param r: arc radius in um, defaults to 10
        :type r: float, optional
        :param n_pts: number of points, defaults to 128
        :type n_pts: int, optional
        :param φ: arc angle in degrees, defaults to 90
        :type φ: float, optional
        :param φ0: arc angle center value in degrees, defaults to 0
        :type φ0: float, optional
        :return: points list
        :rtype: list
        """

        pts = [pya.DPoint(r*np.cos(φ), r*np.sin(φ))
               for φ in np.linspace(np.deg2rad(φ0 - φ/2), np.deg2rad(φ0 + φ/2), n_pts)]

        return pts

    def _produce_body(self):
        """Produces qubit body and stubs 
        """
        dbu = self.layout.dbu

        # Produce circular pad

        body_reg = self._circle_reg_from_radius(self.body_radius)

        # Produce stubs
        stub_end = self.body_radius + 2*self.body_ring_gap + self.coupling_ring_width
        stub_pts = [pya.DPoint(0.0, 0.0),
                    pya.DPoint(0.0, -stub_end)]

        stub_poly = pya.DPath(stub_pts, self.stub_width)
        stub_reg = pya.Region(stub_poly.to_itype(dbu))

        stub_arc_poly = pya.DPath(self._produce_arc_pts(
            r=stub_end, φ=self.stub_angle, φ0=-90), self.stub_width)
        stub_arc_reg = pya.Region(stub_arc_poly.to_itype(dbu))

        stub_reg += stub_arc_reg

        # Produce ring

        # Inner ring
        ring_inner_reg = self._circle_reg_from_radius(
            self.body_radius + self.body_ring_gap)

        # Outer ring
        ring_outer_reg = self._circle_reg_from_radius(
            self.body_radius + self.body_ring_gap + self.coupling_ring_width)

        # stub removal

        stub_cut_path = pya.DPath(
            stub_pts, self.stub_width + 2*self.body_ring_gap)
        stub_cut_reg = pya.Region(stub_cut_path.to_itype(dbu))

        ring_reg = ring_outer_reg - ring_inner_reg - stub_cut_reg

        # Produce circular gap
        body_gap_reg = self._circle_reg_from_radius(
            self.body_radius + self.body_ring_gap + self.coupling_ring_width + self.ring_gnd_gap)

        positive_shape = body_reg + ring_reg + stub_reg

        negative_shape = body_gap_reg

        final_shape = negative_shape - positive_shape

        # Insert shape
        self.cell.shapes(self.get_layer(
            "base_metal_gap_wo_grid")).insert(final_shape)

        # Set probe point
        probe_point = pya.DPoint(0, 0)
        self.refpoints["probe_qb_c"] = probe_point

    def _produce_ring(self):
        """Produces coupling ring
        """

        pass

    def _produce_charge_line(self):
        """Produces charge line
        """
        pass
