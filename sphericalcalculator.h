#ifndef SPHERICALCALCULATOR_H
#define SPHERICALCALCULATOR_H

#include <vector>
#include <string>


const float pi=3.14159265358979323846;
#define SQR(x) (x)*(x)

//limits for maximum L and M calculation. Needed for precalculation of factorials.
const int maxL=20;
const int maxM=20;
const int maxfactorial=maxL+maxM;

struct Point
{
 float x[3];
};

struct Face
{
    float x[3][3];
};

const Point zero={{0,0,0}};

struct Overlay
{
  float intensity;
  int l,m;
};

class SphericalCalculator
{
public:
    SphericalCalculator();
    void setResolution(int _numLayers, int _numAngles);
    void update();

    void saveFile(std::string filename);

    int numAngles;
    int numLayers;

    std::vector< std::vector<Point> > pts; //point array

    std::vector<Face> faces;
    std::vector<Overlay> overlays;

private:
    double factorial [maxfactorial] ;  //factorial[x]=x!
    double nlm [maxL][maxM] ; //normalization constant of the spherical harmonic

    float sphericalHarmonic(int l,int m, float cosphi);

    void calculate();
    void createFaces();
    void addCap( int l);
    void addFace(const Point &a,const Point &b,const Point &c);
};

#endif // SPHERICALCALCULATOR_H
