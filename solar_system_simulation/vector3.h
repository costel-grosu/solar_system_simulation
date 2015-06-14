#ifndef __VECTOR3_H
#define __VECTOR3_H

#include <cmath>
#include <algorithm>

#ifndef INLINE
#ifdef __GNUC__
#define INLINE __attribute__((always_inline))
#else
#define INLINE __forceinline
#endif
#endif

template <typename REAL> 
struct vector3{
	REAL x, y, z;
	INLINE vector3(const vector3 &rhs) : x(rhs.x), y(rhs.y), z(rhs.z) {}
	INLINE vector3(){ x = y = z = REAL(0);}
	INLINE vector3(const REAL &r){x = y = z = r; }
	INLINE vector3(const REAL &_x, const REAL &_y, const REAL &_z){	x = _x;  y = _y;  z = _z;}
	INLINE vector3(const REAL *p){x = p[0]; y = p[1]; z = p[2];}
	INLINE void store(REAL *p){p[0] = x; p[1] = y; p[2] = z;}
	INLINE ~vector3(){}

	INLINE REAL &operator [](int i){ return (&x)[i];}
	INLINE const REAL &operator [](int i) const{ return (&x)[i];}
	template <class real> INLINE operator vector3<real> () const {return vector3<real> (real(x), real(y), real(z));	}
	INLINE operator REAL *(){return &x;}
	INLINE REAL (*toPointer())[3]{return (REAL (*)[3])&x;}
	typedef REAL (*pArrayOfReal3)[3];
	INLINE operator pArrayOfReal3(){return toPointer();}
	INLINE bool are_numbers () const{ return (norm2() >= REAL(0));}
	INLINE REAL norm2() const{return (*this)*(*this);}
	INLINE REAL abs() const{return std::sqrt(norm2());}
	INLINE const vector3<REAL> unit() const{return *this / abs();}

	INLINE const vector3<REAL> operator + (const vector3<REAL> &v) const{return vector3<REAL> (x+v.x, y+v.y, z+v.z);}
	INLINE const vector3<REAL> operator - (const vector3<REAL> &v) const{return vector3<REAL> (x-v.x, y-v.y, z-v.z);}
	INLINE const vector3<REAL> operator * (const REAL &s) const{return vector3<REAL> (x*s, y*s, z*s);}
	INLINE friend const vector3<REAL> operator * (const REAL &s, const vector3<REAL> &v){return v*s;}
	// dot product
	INLINE const REAL operator * (const vector3<REAL> &v) const{return (x*v.x + y*v.y + z*v.z);}
	// vector product
	INLINE const vector3<REAL> operator % (const vector3<REAL> &v) const{return vector3<REAL> (y*v.z - z*v.y, z*v.x - x*v.z, x*v.y - y*v.x);}
	INLINE const vector3<REAL> operator / (const REAL &s) const{REAL r = REAL(1)/s;	return (*this)*r;}
	INLINE const vector3<REAL> &operator = (const vector3<REAL> &v){x = v.x; y=v.y; z=v.z;	return *this;}

	INLINE const vector3<REAL> operator - () const {return vector3<REAL> (-x, -y, -z);}
	INLINE const vector3<REAL> &operator += (const vector3<REAL> &v){*this = *this + v;	return *this;}
	INLINE const vector3<REAL> &operator -= (const vector3<REAL> &v){*this = *this - v;	return *this;}
	INLINE const vector3<REAL> &operator *= (const REAL &s){*this = *this * s;return *this;}
	INLINE const vector3<REAL> &operator /= (const REAL &s){*this = *this / s;	return *this;}

	INLINE friend const vector3<REAL> maxeach (const vector3<REAL> &a, const vector3<REAL> &b){ return vector3<REAL> (std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z));}
	INLINE friend const vector3<REAL> mineach (const vector3<REAL> &a, const vector3<REAL> &b){	return vector3<REAL> (std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z));}
	INLINE const vector3<REAL> abseach(){return vector3<REAL> (std::fabs(x), std::fabs(y), std::fabs(z));}

};

typedef vector3<double> dvec3;
typedef vector3<float>  fvec3;

#undef INLINE

#endif //  __VECTOR3_H
