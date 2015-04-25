# distutils: language = c++
# distutils: sources = cpp/Poisson.cpp
# distutils: extra_compile_args = [-std=c++11]

from libcpp.vector cimport vector
from libcpp cimport bool


cdef extern from "cpp/Poisson.h":
    cdef cppclass sPoint:
        sPoint(float x, float y)
        float x
        float y
        bool IsInRectangle()
        bool IsInCircle()

# cdef sPoint *k = new sPoint(3.0, -2.0)
# print(k.x)
# print(k.y)


cdef extern from "cpp/Poisson.h":
    cdef vector[sPoint] gen_circ_distribution(unsigned seed, float rdisk,
                                              float rscat, float ff,
                                              float distance_border,
                                              bool allow_disks_on_boundary)
    cdef vector[sPoint] gen_rect_distribution(unsigned seed, float w, float h,
                                              float rscat, float ff,
                                              float distance_border,
                                              bool allow_disks_on_boundary)

def generate_circ_poisson_points(seed, rdisk, rscat, ff,
                                 distance_border=0.,
                                 allow_disks_on_boundary=False):
    cdef vector[sPoint] aa

    # TODO use numpy ndarray instead of python list as a return value
    b = []
    aa = gen_circ_distribution(seed, rdisk, rscat, ff, distance_border,
                               allow_disks_on_boundary)
    N=aa.size()
    for i in range(N):
        b.append(aa[i].x + 1j*aa[i].y) # creates the list from the vector
    return b

def generate_rect_poisson_points(seed, w, h, rscat, ff,
                                 distance_border=0.,
                                 allow_disks_on_boundary=False):
    cdef vector[sPoint] aa

    # TODO use numpy ndarray instead of python list as a return value
    b = []
    aa = gen_rect_distribution(seed, w, h, rscat, ff, distance_border,
                               allow_disks_on_boundary)
    N=aa.size()
    for i in range(N):
        b.append(aa[i].x + 1j*aa[i].y) # creates the list from the vector
    return b
