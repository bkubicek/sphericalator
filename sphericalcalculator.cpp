#include "sphericalcalculator.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>


SphericalCalculator::SphericalCalculator()
{
    setResolution(130,130);


     //prepare factorial and nlm
     factorial[0]=1;
     factorial[1]=1;
     for(int i=2;i<maxfactorial;i++)
        factorial[i]=factorial[i-1]*i;

     for(int l=0;l<maxL;l++)
     for(int m=0;m<maxM;m++)
        nlm[l][m]=sqrt((2*l+1)*factorial[l-m]/factorial[l+m]);
}

void SphericalCalculator::setResolution(int _numLayers, int _numAngles)
{
    numLayers=_numLayers;
    numAngles=_numAngles;
    pts.resize( numLayers );
    for(int l=0;l<numLayers;l++)
      pts[l].resize(numAngles);
    faces.resize(0);
}



using namespace std;


#define GSL_DBL_EPSILON        2.2204460492503131e-16
#define GSL_DBL_MIN        2.2250738585072014e-308
#define GSL_LOG_DBL_MIN   -7.0839641853226408e+02

double gsl_sf_legendre_Plm(const int l,const int m, const double x)
{
  /* If l is large and m is large, then we have to worry
   * about overflow. Calculate an approximate exponent which
   * measures the normalization of this thing.
   */

    double dif = l-m;
      double sum = l+m;
    double exp_check = 0.5 * log(2.0*l+1.0)
                       + 0.5 * dif * (log(dif)-1.0)
                       - 0.5 * sum * (log(sum)-1.0);


  if(m < 0 || l < m || x < -1.0 || x > 1.0) {
    //DOMAIN_ERROR(result);
      return 0;
  }
  else if(exp_check < GSL_LOG_DBL_MIN + 10.0){

    //OVERFLOW_ERROR(result);
  return 0;
  }
  else {
    /* Account for the error due to the
     * representation of 1-x.
     */


    double p_mm;     /* P_m^m(x) */
    double p_mmp1;   /* P_{m+1}^m(x) */

    /* Calculate P_m^m from the analytic result:
     *          P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^(m/2) , m > 0
     */
    p_mm = 1.0;
    if(m > 0){
      double root_factor = sqrt(1.0-x)*sqrt(1.0+x);
      double fact_coeff = 1.0;
      int i;
      for(i=1; i<=m; i++) {
        p_mm *= -fact_coeff * root_factor;
        fact_coeff += 2.0;
      }
    }

    /* Calculate P_{m+1}^m. */
    p_mmp1 = x * (2*m + 1) * p_mm;

    if(l == m){
      return p_mm;


    }
    else if(l == m + 1) {
      return   p_mmp1;
    }
    else{
      double p_ell = 0.0;
      int ell;

      /* Compute P_l^m, l > m+1 by upward recurrence on l. */
      for(ell=m+2; ell <= l; ell++){
        p_ell = (x*(2*ell-1)*p_mmp1 - (ell+m-1)*p_mm) / (ell-m);
        p_mm = p_mmp1;
        p_mmp1 = p_ell;
      }

      return p_ell;

    }
  }
  return 0;
}


float SphericalCalculator::sphericalHarmonic(int l,int m, float cosphi)
{
 return nlm[l][m]*gsl_sf_legendre_Plm (l, m, cosphi);
}



void SphericalCalculator::calculate()
{
   #pragma omp parallel for
   for(int l=0;l<numLayers;l++)
   {
     //  double sph[maxL][maxM];
       std::vector<double> sphf;
       sphf.resize(overlays.size());
     float cosphi=cos(pi*l/float(numLayers-1));
     float sinphi=sin(pi*l/float(numLayers-1));
     for(int ii=0;ii<overlays.size();ii++)
         sphf[ii]=sphericalHarmonic(overlays[ii].l,overlays[ii].m,cosphi);
     /*
     for(int ll=0;ll<maxL;ll++)
     for(int m=0;m<=ll;m++)
       sph[ll][m]=sphericalHarmonic(ll,m,cosphi);
    */


     for(int i=0;i<numAngles;i++)
     {
        Point &p=pts[l][i];
       //float h=l*height/float(numLayer)-R;

       float alpha=2*i/float(numAngles)*pi;


       float r=0;
       for(int i=0;i<overlays.size();i++)
       {
           Overlay &o=overlays[i];
           //r+=o.intensity*sph[o.l][o.m]*cos(o.m*alpha);
           r+=o.intensity*sphf[i]*cos(o.m*alpha);
       }
         //0.1*sph[5][5]*cos(5*alpha)+
         //0.1*sph[5][2]*cos(2*alpha);
        // 0.25*sph[19][n]*cos(n*alpha);
//          rr*0.05*exp(-SQR((h)/3))*sin(alpha*3)+
// 	 rr*0.05*exp(-SQR((h)/3))*sin(alpha*8)+
// 	 rr*0.05*exp(-SQR((h-R/2.)/2))*sin(alpha*16)+
// 	 rr*0.05*exp(-SQR((h+R/2.)/2))*cos(alpha*16);


       if(r<0)r=0;
       p.x[0]=r*cos(alpha)*sinphi;
       p.x[1]=r*sin(alpha)*sinphi;
       p.x[2]=r*cosphi;

     }
   }
}

void putpoint(ostream &out, const Point &p)
{
  out<<"vertex ";
  for(int j=0;j<3;j++)
    out<<p.x[j]<<" ";
  out<<endl;
}

void SphericalCalculator::addCap( int l)
{
  bool bottom=false;
  if(l==0) bottom=true;

  Point cpos={{0,0,0}};
  for(int i=0;i<numAngles;i++)
  {
      cpos.x[0]+=pts[l][i].x[0];
      cpos.x[1]+=pts[l][i].x[1];
      cpos.x[2]+=pts[l][i].x[2];
  }
  cpos.x[0]*=1/float(numAngles);
  cpos.x[1]*=1/float(numAngles);
  cpos.x[2]*=1/float(numAngles);

  for(int i=0;i<numAngles-1;i++)
  {
      if(bottom)
        addFace(cpos,pts[l][i],pts[l][i+1]);
      else
        addFace(cpos,pts[l][i+1],pts[l][i]);
  }
  if(bottom)
    addFace(cpos,pts[l][numAngles-1],pts[l][0]);
  else
    addFace(cpos,pts[l][0],pts[l][numAngles-1]);
}

void SphericalCalculator::addFace(const Point &a,const Point &b,const Point &c)
{
    Face f;
    f.x[0][0]=a.x[0];
    f.x[0][1]=a.x[1];
    f.x[0][2]=a.x[2];

    f.x[1][0]=b.x[0];
    f.x[1][1]=b.x[1];
    f.x[1][2]=b.x[2];

    f.x[2][0]=c.x[0];
    f.x[2][1]=c.x[1];
    f.x[2][2]=c.x[2];

    faces.push_back(f);
}

void SphericalCalculator::createFaces()
{
    faces.clear();
    for(int l=0;l<numLayers-1;l++)
    {
      for(int i=0;i<numAngles-1;i++)
      {
        addFace(pts[l][i], pts[l+1][i],   pts[l+1][i+1]);
        addFace(pts[l][i], pts[l+1][i+1], pts[l][i+1]  );
      }
      addFace(pts[l][0], pts[l+1][numAngles-1],   pts[l+1][0]);
      addFace(pts[l][0], pts[l][numAngles-1],     pts[l+1][numAngles-1]);


    }
    addCap(0);
    addCap(numLayers-1);
}

void SphericalCalculator::saveFile(std::string  filename)
{
  cerr<<"exporting"<<endl;

  fstream out(filename.c_str(),fstream::out);
  out<<"solid thing"<<endl;
  for(int i=0;i<(int)faces.size();i++)
  {
    out<<"faced normal 0 0 0\nouter loop\n";
    out<<"vertex ";
    out<<faces[i].x[0][0]<<" ";
    out<<faces[i].x[0][1]<<" ";
    out<<faces[i].x[0][2]<<endl;
    out<<"vertex ";
    out<<faces[i].x[1][0]<<" ";
    out<<faces[i].x[1][1]<<" ";
    out<<faces[i].x[1][2]<<endl;
    out<<"vertex ";
    out<<faces[i].x[2][0]<<" ";
    out<<faces[i].x[2][1]<<" ";
    out<<faces[i].x[2][2]<<endl;
    out<<"endloop\nendfacet\n";
  }
  out<<"endsolid thing"<<endl;
  out.close();
}

void SphericalCalculator::update()
{
    calculate();
    createFaces();
}
