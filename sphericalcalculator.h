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
    float n[3];
};

struct Edge
{
  Point p[2];
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
    void flatbottom();
    float overhangstart;

    void saveFile(std::string filename);

    int numAngles;
    int numLayers;
    float minz;
    float radius;
    float supportgap,mygap; //mygap=support gap brokne down to r=1;
    bool support;
    float flatheight;

    std::vector< std::vector<Point> > pts; //point array

    std::vector<Face> faces;
    std::vector<Face> overhangs;
    std::vector<Overlay> overlays;
    
    std::vector<Edge> edges;
    std::vector<int> count;
    void addEdge(const Point &a,const Point &b);

private:
    double factorial [maxfactorial] ;  //factorial[x]=x!
    double nlm [maxL][maxM] ; //normalization constant of the spherical harmonic

    float sphericalHarmonic(int l,int m, float cosphi);

    void calculate();
    void createFaces();
    void addCap( int l);
    void addFace(const Point &a,const Point &b,const Point &c);
    void addFace(std::vector<Face> &faces,const Point &a,const Point &b,const Point &c);
    void createSupport();
};

#endif // SPHERICALCALCULATOR_H
