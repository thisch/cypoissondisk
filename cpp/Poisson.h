#include <vector>


struct sPoint {
  sPoint() : x(0), y(0), m_Valid(false) {}
  sPoint(float X, float Y) : x(X), y(Y), m_Valid(true) {}
  float x;
  float y;
  bool m_Valid;
  //
  bool IsInRectangle() const;
  //
  float R() const {
    float fx = 2*(x - 0.5f);
    float fy = 2*(y - 0.5f);
    return sqrt(fx*fx + fy*fy);
  }
  //
 bool IsInCircle() const;
};


std::vector<sPoint> gen_circ_distribution(unsigned, float, float, float, float, bool);
std::vector<sPoint> gen_rect_distribution(unsigned, float, float, float, float, float, bool);
