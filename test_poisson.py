import logging
from unittest import TestCase

import pytest
import numpy as np
from scipy.spatial import distance

from poissondisk import generate_circ_poisson_points as gpp
from poissondisk import generate_rect_poisson_points as grpp

LG = logging.getLogger('cypoisson')


class TestBase(TestCase):

    def highlight(self, text):
        LG.info(text.center(80, '#'))

    def _plot(self, rect=False):
        points = np.array(self.points)
        lp = len(points)
        import matplotlib
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle, Rectangle
        from matplotlib.collections import PatchCollection

        fig = plt.figure()
        ax = fig.add_subplot(111)

        patches = []
        for p in points:
            circle = Circle((p.real, p.imag), self.rscat)
            patches.append(circle)

        try:
            w = h = 2*self.rdisk
        except AttributeError:
            w, h = self.w, self.h
        if rect:
            patches.append(Rectangle((-w/2., -h/2.), w, h))
        else:
            patches.append(Circle((0, 0), self.rdisk))

        p = PatchCollection(patches, cmap=matplotlib.cm.rainbow,
                            alpha=0.4)
        ax.add_collection(p)

        pad = 0.05*w
        ax.set_xlim(-w/2. - pad, w/2. + pad)
        ax.set_ylim(-h/2. - pad, h/2. + pad)
        ax.grid(True)
        ax.set_xlabel('x [um]')
        ax.set_ylabel('y [um]')
        ax.set_aspect('equal')

        if rect:
            ax.set_title('%d scatterers (w %g um, h %g um, ff=%g)' % (
                lp, w, h, self.ff))
        else:
            ax.set_title('%d scatterers (diameter %g um, ff=%g)' % (
                lp, self.rdisk*2, self.ff))
        plt.show()


class TestPoisson(TestBase):
    """
    TODO
    unterschiedliche seeds muessen unterschiedliche konfigs erzeugen
    gleiche seeds gleiche konfigs
    distance_border
    max distance ~ ddisk-distanceborder

    """

    def test_1(self):
        self.rdisk = 350.
        self.rscat = 10.
        desired_ff = 0.4
        self.points = gpp(0, rdisk=self.rdisk, rscat=self.rscat,
                          ff=desired_ff)
        lp = len(self.points)
        self.assertGreater(lp, 0)
        ff = lp*self.rscat**2./self.rdisk**2
        LG.info('desired ff %g, actual ff %g', desired_ff, ff)
        self.assertLess(abs(desired_ff - ff), 0.03)
        self.ff, self.desired_ff = ff, desired_ff

    def test_1_rectangle(self):
        self.h = self.w = 700.
        self.rscat = 10.
        desired_ff = 0.45
        self.points = grpp(1, w=self.w, h=self.h, rscat=self.rscat,
                           ff=desired_ff)
        lp = len(self.points)
        self.assertGreater(lp, 0)
        ff = lp*np.pi*self.rscat**2./(self.w*self.h)
        LG.info('desired ff %g, actual ff %g', desired_ff, ff)
        self.assertLess(abs(desired_ff - ff), 0.05)
        self.ff, self.desired_ff = ff, desired_ff

    @pytest.mark.skipif('not config.getvalue("interactive")')
    def test_1_plot(self):
        self.test_1()
        self._plot()

    @pytest.mark.skipif('not config.getvalue("interactive")')
    def test_1_rectangle_plot(self):
        self.test_1_rectangle()
        self._plot(rect=True)

    def test_2(self):
        self.rdisk = 350.
        self.rscat = 10.
        prevlp = 0
        for desired_ff in np.linspace(0.05, 0.4, 8):
            points = gpp(1, rdisk=self.rdisk, rscat=self.rscat, ff=desired_ff)
            ff = len(points)*self.rscat**2./self.rdisk**2
            LG.info('desired ff %g, actual ff %g', desired_ff, ff)
            self.assertLess(abs(desired_ff - ff), 0.03)
            self.assertGreater(len(points), prevlp)
            levlp = len(points)

    def test_3(self):
        """
        collision check
        check that all points lie within a the disk
        """
        self.rdisk = 7.5
        self.rscat = 0.1
        desired_ff = 0.30
        self.points = np.array(gpp(1, rdisk=self.rdisk, rscat=self.rscat,
                                   ff=desired_ff))
        self.assertGreater(len(self.points), 0)
        ff = len(self.points)*self.rscat**2./self.rdisk**2
        LG.info('desired ff %g, actual ff %g', desired_ff, ff)
        self.assertLess(abs(desired_ff - ff), 0.03)
        self.ff, self.desired_ff = ff, desired_ff

        assert all(np.abs(self.points) < self.rdisk - self.rscat)

        # self._plot()

        p = np.c_[self.points.real, self.points.imag]
        dists = distance.cdist(p, p, 'euclidean')
        mindist = dists[dists > 0].min()
        maxdist = dists.max()
        self.assertGreater(mindist, 2*self.rscat)

        LG.info("mindist: %g\tnet mindist: %g\tmaxdist: %g",
                mindist, mindist - 2*self.rscat, maxdist)
        self.assertGreater(maxdist, self.rdisk*2*0.9)

    def test_3_loop(self):
        """
        collision check
        check that all points lie within a the disk
        """
        self.rdisk = 7.5  # ddisk = 15um
        self.rscat = 0.05  # dscat = 100nm
        for desired_ff in [0.05, 0.10, 0.15, 0.2, 0.25, 0.3, 0.35]:
            self.highlight("DESIRED FF %g" % desired_ff)
            self.points = np.array(gpp(1, rdisk=self.rdisk, rscat=self.rscat,
                                       ff=desired_ff))
            self.assertGreater(len(self.points), 0)
            ff = len(self.points)*self.rscat**2./self.rdisk**2
            LG.info('desired ff %g, actual ff %g', desired_ff, ff)
            self.assertLess(abs(desired_ff - ff), 0.03)
            self.ff, self.desired_ff = ff, desired_ff

            assert all(np.abs(self.points) < self.rdisk - self.rscat)

            # self._plot()
            p = np.c_[self.points.real, self.points.imag]
            dists = distance.cdist(p, p, 'euclidean')
            mindist = dists[dists > 0].min()
            maxdist = dists.max()
            self.assertGreater(mindist, 2*self.rscat)

            LG.info("FF %g:\tmindist: %g\tnet mindist: %g\tmaxdist: %g",
                    self.ff, mindist, mindist - 2*self.rscat, maxdist)
            self.assertGreater(maxdist, self.rdisk*2*0.9)

    def test_4(self):
        """
        test distributions where disks intersecting the boundary are allowed
        """
        self.rdisk = 350.
        self.rscat = 10.
        desired_ff = 0.4
        self.points = gpp(0, rdisk=self.rdisk, rscat=self.rscat,
                          ff=desired_ff, allow_disks_on_boundary=True)
        lp = len(self.points)
        self.assertGreater(lp, 0)
        ff = lp*self.rscat**2./self.rdisk**2  # TODO this is not the correct
                                              # expression for the filling
                                              # fraction as some of the
                                              # scatterers lie outside the
                                              # disk
        LG.info('desired ff %g, actual ff %g', desired_ff, ff)
        self.assertLess(abs(desired_ff - ff), 0.01)
        self.ff, self.desired_ff = ff, desired_ff

        expr = [abs(p) > self.rdisk-self.rscat for p in self.points]
        LG.info("%d scats lie on the boundary", sum(expr))
        assert any(expr)

    @pytest.mark.skipif('not config.getvalue("interactive")')
    def test_4_plot(self):
        self.test_4()
        self._plot()
