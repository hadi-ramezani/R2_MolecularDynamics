#ifndef VECTOR_H
#define	VECTOR_H

#include <iostream>
#include <math.h>
#include <stdio.h>


class Vector {
   public:
     double x,y,z;
     
     Vector(void) {         // default is to create a 0 vector
       x = y = z = 0.0;
     }
     Vector( const Vector &v2) { // Vector x = another_vector
       x = v2.x;
       y = v2.y;
       z = v2.z;
     }
     Vector( double newx, double newy, double newz) {
       x = newx;
       y = newy;
       z = newz;
     }
     ~Vector( void) { }

     //  v1 = v2;
     Vector& operator=(const Vector &v2) {
       x = v2.x;
       y = v2.y;
       z = v2.z;
       return *this;
     }

     //  v1 = const;
     Vector& operator=(const double &v2) {
       x = v2;
       y = v2;
       z = v2;
       return *this;
     }

     //  v1 += v2;
     void operator+=(const Vector &v2) {
       x += v2.x;
       y += v2.y;
       z += v2.z;
     }

     // v1 -= v2;
     void operator-=(const Vector &v2) {
       x -= v2.x;
       y -= v2.y;
       z -= v2.z;
     }

     // v1 *= const
     void operator*=(const double &v2) {
       x *= v2;
       y *= v2;
       z *= v2;
     }

     // v1 /= const
     void operator/=(const double& v2) {
       x /= v2;
       y /= v2;
       z /= v2;
     }

     friend int operator == (const Vector& v1, const Vector& v2) {
       return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
     }
     friend int operator != (const Vector& v1, const Vector& v2) {
       return v1.x != v2.x || v1.y != v2.y || v1.z != v2.z;
     }

     // addition of two vectors
     friend Vector operator+(const Vector& v1, const Vector& v2) {
       return Vector( v1.x+v2.x, v1.y+v2.y, v1.z+v2.z);
     }

     // negation
     friend Vector operator-(const Vector &v1) {
       return Vector( -v1.x, -v1.y, -v1.z);
     }

     // subtraction
     friend Vector operator-(const Vector &v1, const Vector &v2) {
       return Vector( v1.x-v2.x, v1.y-v2.y, v1.z-v2.z);
     }
     // inner ("dot") product
     friend double operator*(const Vector &v1, const Vector &v2) {
       return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
     }
     // scalar product
     friend Vector operator*(const double &f, const Vector &v1) {
       return Vector(f*v1.x, f*v1.y, f*v1.z);
     }
     // scalar product
     friend Vector operator*(const Vector &v1, const double &f) {
       return Vector(f*v1.x, f*v1.y, f*v1.z);
     }
     // division by a scalar
     friend Vector operator/(const Vector &v1, const double &f) {
       return Vector(v1.x/f, v1.y/f, v1.z/f);
     }
     
     // return the norm
     inline double length(void) const {
       return sqrt(x*x+y*y+z*z);
     }
     
     double length2(void) const {
       return (x*x + y*y + z*z);
     }

     // return the unit vector in the same direction
     Vector unit(void) const {
       return Vector(x, y, z)/length();
     }
     
     // print out
     friend std::ostream& operator<<(std::ostream& strm, const Vector &v1) {
       strm << "( "<< v1.x << ", " << v1.y << ", " << v1.z << ')';
       return strm;
     }

     // one cross product  self = cross(v1, v2)
     void cross(const Vector &v1, const Vector &v2) {
       x= v1.y*v2.z-v2.y*v1.z;
       y= v2.x*v1.z-v1.x*v2.z;
       z= v1.x*v2.y-v2.x*v1.y;
     }

     // multiplying a cross product by a scalar is very common
     // one cross product  v3 = k*cross(v1, v2)
     inline friend Vector cross(const double k, const Vector &v1, const Vector &v2) {
       return Vector( k*(v1.y*v2.z-v2.y*v1.z),
                      // k*(-v1.x*v2.z+v2.x*v1.z),
                      k*(v2.x*v1.z-v1.x*v2.z),
                      k*(v1.x*v2.y-v2.x*v1.y) );
     }

     // add a vector to this vector
     void add(const Vector &v1) {
       x+=v1.x; y+=v1.y; z+=v1.z;
     }

     // subtract the vector from this one
     void sub(const Vector &v1) {
       x-=v1.x; y-=v1.y; z-=v1.z;
     }

     // add a constant factor to each element of a vector
     void add_const(double c)
     {
	x+=c;
	y+=c;
	z+=c;
     }

     // rescale everything by a scalar -- V = a*V
     void mult(double f) {
       x*=f; y*=f; z*=f;
     }

     // divide each element by a scalar
     void div(double f) {
       x/=f; y/=f; z/=f;
     }

     // returns (*this) * V2
     double dot(const Vector &v2) {
       return x*v2.x + y*v2.y + z*v2.z;
     }

     // set the vector based on a string.  If bad, return 0
     // the string can be in the form "x y z" or "x, y, z"
     int set(const char *s) {
	double a[3];    // read into doubles, since I don't know what
	char tmp[100];  // a "double" is in real life
	// cheap way to get commas, etc.  a poor regex
       int i=sscanf(s, "%lf%99[ \t,]%lf%99[ \t,]%lf%99s",
                    a, tmp, a+1, tmp, a+2, tmp);
       if (i != 5) return 0;
       const char *t = s;       // now count commas (for "1,,,,2,  , 3")
       int flg = 0;                 // and check for "1 2,,3"
       i = 0;
       for (;*t;t++) {
          if (*t == ',') { 
             if (flg == 0) {   // expecting non-whitespace
                return 0;  //    so error
             }
             flg = 0;          // now expect non-whitespace
             i++;              // and increment comma counter
          }
          else if (*t != ' ' && *t != '\t') {  // got non-whitespace
             flg = 1;          // so the next can be whitespace or commas
          }
       }
       if (i == 0 || i == 2) {  // allow "1 2 3" or "1, 2,3" forms
          x = a[0]; y = a[1]; z = a[2];
          return 1;
       }
       return 0;
     }
};

#endif	/* VECTOR_H */


