/**
 * \file Poisson.cpp
 * \brief
 *
 * Poisson Disk Points Generator
 *
 * \version 1.1.1
 * \date 23/05/2014
 * \author Sergey Kosarevsky, 2014
 * \author support@linderdaum.com   http://www.linderdaum.com
 *http://blog.linderdaum.com
 */

// Fast Poisson Disk Sampling in Arbitrary Dimensions
// http://people.cs.ubc.ca/~rbridson/docs/bridson-siggraph07-poissondisk.pdf

// Implementation based on
// http://devmag.org.za/2009/05/03/poisson-disk-sampling/

#include <iostream>
#include <vector>
#include <random>
#include <stdint.h>
#include <time.h>
#include <fstream>
#include <memory.h>

#include "Poisson.h"

///////////////// User selectable parameters ///////////////////////////////

double rscat = 10.;
double ff = 0.35;
double rdisk = 350;
double distance_border;
double w;
double h;
int NumPoints;
float MinDistance;

int k = 30; // refer to bridson-siggraph07-poissondisk.pdf for details

////////////////////////////////////////////////////////////////////////////

std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> dis(0.0, 1.0);


bool sPoint::IsInCircle() const {
    float fx = 2*(x - 0.5f);
    float fy = 2*(y - 0.5f);
    return sqrt(fx * fx + fy * fy) < (1-(rscat+distance_border)/rdisk);
}

bool sPoint::IsInRectangle() const {
    return
        x >= (rscat+distance_border)/w &&
        y >= (rscat+distance_border)/h &&
        x <= (1-(rscat+distance_border)/w) &&
        y <= (1-(rscat+distance_border)/h);
        // x >= 0 && y >= 0 && x <= 1 && y <= 1;
}


struct sGridPoint {
  sGridPoint(int X, int Y) : x(X), y(Y) {}
  int x;
  int y;
};

float RandomFloat() { return static_cast<float>(dis(gen)); }

float GetDistance(const sPoint &P1, const sPoint &P2) {
  return sqrt((P1.x - P2.x) * (P1.x - P2.x) + (P1.y - P2.y) * (P1.y - P2.y));
}

sGridPoint ImageToGrid(const sPoint &P, float CellSize) {
  return sGridPoint((int)(P.x / CellSize), (int)(P.y / CellSize));
}

struct sGrid {
  sGrid(int W, int H, float CellSize) : m_W(W), m_H(H), m_CellSize(CellSize) {
    m_Grid.resize(m_H);

    for (auto i = m_Grid.begin(); i != m_Grid.end(); i++) {
      i->resize(m_W);
    }
  }
  void Insert(const sPoint &P) {
    sGridPoint G = ImageToGrid(P, m_CellSize);
    m_Grid[G.x][G.y] = P;
  }
  bool IsInNeighbourhood(sPoint Point, float MinDist, float CellSize) {
    sGridPoint G = ImageToGrid(Point, CellSize);

    // number of adjucent cells to look for neighbour points
    const int D = 5;

    // scan the neighbourhood of the point in the grid
    for (int i = G.x - D; i < G.x + D; i++) {
      for (int j = G.y - D; j < G.y + D; j++) {
        if (i >= 0 && i < m_W && j >= 0 && j < m_H) {
          sPoint P = m_Grid[i][j];

          if (P.m_Valid && GetDistance(P, Point) < MinDist) {
            return true;
          }
        }
      }
    }

    return false;
  }

private:
  int m_W;
  int m_H;
  float m_CellSize;

  std::vector<std::vector<sPoint>> m_Grid;
};

sPoint PopRandom(std::vector<sPoint> &Points) {
  std::uniform_int_distribution<> dis(0, Points.size() - 1);
  int Idx = dis(gen);
  sPoint P = Points[Idx];
  Points.erase(Points.begin() + Idx);
  return P;
}

sPoint GenerateRandomPointAround(const sPoint &P, float MinDist) {
  // start with non-uniform distribution
  float R1 = RandomFloat();
  float R2 = RandomFloat();

  // radius should be between MinDist and 2 * MinDist
  float Radius = MinDist * (R1 + 1.0f);

  // random angle
  float Angle = 2 * M_PI * R2;

  // the new point is generated around the point (x, y)
  float X = P.x + Radius * cos(Angle);
  float Y = P.y + Radius * sin(Angle);

  return sPoint(X, Y);
}

std::vector<sPoint> GeneratePoissonPoints(float MinDist, int NewPointsCount,
                                          size_t NumPoints, bool Circle,
                                          bool allow_disks_on_boundary=false) {
  std::vector<sPoint> SamplePoints;
  std::vector<sPoint> ProcessList;

  // create the grid
  float CellSize = MinDist / sqrt(2.0f);

  int GridW = (int)ceil(1.0f / CellSize);
  int GridH = (int)ceil(1.0f / CellSize);

  sGrid Grid(GridW, GridH, CellSize);

  if (allow_disks_on_boundary) {
      rscat = 0.;
  }

  sPoint FirstPoint = sPoint(RandomFloat(), RandomFloat());
  do {
    FirstPoint = sPoint(RandomFloat(), RandomFloat());
  } while(Circle ? !FirstPoint.IsInCircle() : !FirstPoint.IsInRectangle());
  // std::cout << "FP " << FirstPoint.R()
  //           << " " << FirstPoint.IsInCircle()
  //           << " x " << FirstPoint.x
  //           << " y " << FirstPoint.y << std::endl;

  // update containers
  ProcessList.push_back(FirstPoint);
  SamplePoints.push_back(FirstPoint);
  Grid.Insert(FirstPoint);

  // generate new points for each point in the queue
  while (!ProcessList.empty() && SamplePoints.size() < NumPoints) {
    // a progress indicator, kind of
    if (SamplePoints.size() % 100 == 0)
      std::cout << ".";

    sPoint Point = PopRandom(ProcessList);

    for (int i = 0; i < NewPointsCount; i++) {
      sPoint NewPoint = GenerateRandomPointAround(Point, MinDist);

      bool Fits = Circle ? NewPoint.IsInCircle() : NewPoint.IsInRectangle();

      if (Fits && !Grid.IsInNeighbourhood(NewPoint, MinDist, CellSize)) {
        ProcessList.push_back(NewPoint);
        SamplePoints.push_back(NewPoint);
        Grid.Insert(NewPoint);
        continue;
      }
    }
  }

  std::cout << std::endl << std::endl;

  return SamplePoints;
}

void PrintBanner() {
  std::cout << "Poisson disk points generator" << std::endl;
  std::cout << "Sergey Kosarevsky, 2014" << std::endl;
  std::cout << std::endl;

  std::cout << "Numpoints " << NumPoints << std::endl;
  std::cout << "dscat " << rscat * 2 << std::endl;
  std::cout << "MinDist " << MinDistance * rdisk * 2 << std::endl;
}

std::vector<sPoint> gen_circ_distribution(unsigned seed, float rdisk_,
                                          float rscat_, float ff_, float distance_border_,
                                          bool allow_disks_on_boundary=false) {
  gen.seed(seed);
  rdisk = rdisk_;
  rscat = rscat_;
  ff = ff_;
  distance_border = distance_border_;

  // TODO do I need to take distance_border into account here?
  NumPoints =
      static_cast<int>(rdisk * rdisk * ff /
                       (rscat * rscat)); // minimal number of points to generate
  MinDistance = 1. / sqrt(float(2 * NumPoints));

  bool circle = true; // all points lie within a circle

  PrintBanner();

  std::vector<sPoint> Points =
      GeneratePoissonPoints(MinDistance, k, NumPoints, circle, allow_disks_on_boundary);
  std::cout <<Points.size() << " Points generated" << std::endl;

  double netrdisk = rdisk - distance_border;
  for (auto& p : Points) {
      p.x = p.x * 2 * netrdisk - netrdisk;
      p.y = p.y * 2 * netrdisk - netrdisk;
  }
  // collision detection has to be done on the client side

  return Points;
}

std::vector<sPoint> gen_rect_distribution(unsigned seed, float w_, float h_,
                                          float rscat_, float ff_, float distance_border_,
                                          bool allow_disks_on_boundary) {
  gen.seed(seed);
  rscat = rscat_;
  ff = ff_;
  distance_border = distance_border_;
  w = w_;
  h = h_;

  // TODO do I need to take distance_border into account here?
  NumPoints =
      static_cast<int>(w * h * ff /
                       (rscat * rscat * M_PI)); // minimal number of points to generate
  MinDistance = 1. / sqrt(float(1.55*NumPoints));

  bool circle = false; // all points lie within a rectangle

  PrintBanner();

  std::vector<sPoint> Points =
      GeneratePoissonPoints(MinDistance, k, NumPoints, circle, allow_disks_on_boundary);
  std::cout <<Points.size() << " Points generated" << std::endl;

  double netw = w - distance_border;
  double neth = h - distance_border;
  for (auto& p : Points) {
      p.x = p.x * netw - netw/2.;
      p.y = p.y * neth - neth/2.;
  }
  // collision detection has to be done on the client side

  return Points;
}
