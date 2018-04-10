/**
* Camera.h
* (c)2008 Takashi Michikawa
*/
#ifndef __CAMERA_H__
#define __CAMERA_H__ 1
#include <vecmath.h>
#include <vector>
#include <algorithm>
#include <iostream>
using std::cerr;
using std::operator<<;
using std::endl;

using namespace kh_vecmath;

class Camera
{
private:

	Vector3d _init_center;		
	double	 _init_translate;   
	double   _init_radius;		
	double   _init_fov;
	
	Vector3d _center;		
	double	 _translate;    
	double   _radius;		
	double   _fov;	

	Quat4d   _q;            ///< ‰ñ“]

public:
	Camera(const double radius = 1.0, const Vector3d& center = Vector3d(0.0, 0.0, 0.0), const Quat4d& rotation = Quat4d(0.0, 0.0, 0.0, 1.0), const double fov = 45.0 )
	{
		setFOV(fov);
		setRadius(radius);
		setCenter(center);
		setQuaternion(rotation);
		_translate = this->get_translate(fov,radius);
	};
	Camera(const Camera& d)
	{
		copy(d);
	}

	Camera& operator = ( const Camera& d)
	{
		copy(d);
		return *this;
	}

	void setInit(void)
	{
		_init_radius = _radius;
		_init_center = _center;
		_init_fov    = _fov;
		_init_translate = _translate;
	}

	void init()
	{
		_radius = _init_radius;
		_center = _init_center;
		_translate = _init_translate;
		_fov = _init_fov;
		_q.set(0.0,0.0,0.0,1.0);
	}

	void applyTransform()
	{
		Matrix4d r; // rotation matrix
		r.setIdentity();
		r.set(_q);

		Point3d eye(0.0, 0.0, _translate);
		Point3d up(0.0, 1.0, 0.0);
		r.transform(&eye);
		r.transform(&up);
		eye.add(_center);

		::gluLookAt(eye.x, eye.y, eye.z, _center.x, _center.y, _center.z, up.x, up.y, up.z);
	}


	void
		updatePerspective()
	{
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();             
		int viewport[4];
		glGetIntegerv(GL_VIEWPORT,viewport);

		double zNear = _translate - _init_radius;
		double zFar  = _translate + _init_radius;
		if(zNear < 0 ) zNear = 0.01;
		gluPerspective(45, viewport[2] * 1.0 / viewport[3],  zNear, zFar);	
	}

	void zoom(const double d)
	{
		if(_translate + d < 0.1) _translate = 0.1;
		else _translate += d;

		_radius = _translate *  sin (_fov / 360.0 * 3.14159) ;
	}

	void rotate(const int oldx, const int oldy, const int newx, const int newy)
	{
		int viewport[4];
		glGetIntegerv(GL_VIEWPORT,viewport);
		const int W = viewport[2];
		const int H = viewport[3];
		const int R = std::min(W,H);

		Quat4d lastq = this->get_rotation(
			(2.0 * oldx - W) / R,    (H - 2.0 * oldy) / R,
			(2.0 * newx - W) / R,    (H - 2.0 * newy) / R
			);
		_q.mul(lastq);
	}
	void move(const double dradius, const Vector3d& dcenter, Quat4d& dq, const double dfov)
	{
		_fov    += dfov;
		_radius += dradius;
		_center += dcenter;
		_translate = get_translate(_fov, _radius);
		_q.mul(dq, _q);
	}

	void setFOV(const double fov)
	{
		if(fov < 10.0)
		{
			cerr<<"invalid fov : "<<fov<<". set to 10"<<endl;
			_fov = 10.0;
		}
		if(fov >= 150.0)
		{
			cerr<<"invalid fov : "<<fov<<". set to 150"<<endl;
			_fov = 150.0;
		}
		_fov = fov;
	};

	double getFOV(void) const
	{
		return _fov;
	}

	void setRadius(const double radius)
	{

		_radius = radius;

	}
	double getRadius(void) const
	{
		return _radius;
	}

	void setCenter(const Vector3d& center)
	{
		_center = center;
	}
	Vector3d getCenter(void) const
	{
		return _center;
	}


	Quat4d getQuaternion(void) const
	{
		return _q;
	}
	void  setQuaternion(const Quat4d& dq) 
	{
		_q = dq;
	}

	Vector3d getRay()const
	{
		Matrix4d r; // rotation matrix
		r.setIdentity();
		r.set(_q);

		Vector3d newRay(0.0, 0.0, 1.0);
		r.transform(&newRay);
		return newRay;
	}

private:

	void project_to_sphere(const double& radius, Vector3d& p) const
	{
		p.z = 0;
		const double d = p.lengthSquared(); //p.x^2 + p.y^2
		const double r = radius * radius;
		if (d < r)	p.z = sqrt(r - d);
		else		p.scale(radius / p.length());
		return;
	}
	Quat4d
		get_rotation(const double& oldx, const double& oldy, const double& newx, const double& newy, const double radius = 0.8) const
	{

		Vector3d oldp(oldx, oldy, 0);
		Vector3d newp(newx, newy, 0);

		if (oldp.epsilonEquals(newp, 1.0e-16 ) )
		{
			return Quat4d(0,0,0,1);
		}

		project_to_sphere(radius, oldp);
		project_to_sphere(radius, newp);
		Vector3d axis;
		axis.cross(newp,oldp);
		const double theta= oldp.angle(newp);
		Quat4d   retquat;
		retquat.set(AxisAngle4d(axis, theta));
		return retquat;
	}

	double get_translate(const double fov, const double radius) const
	{
		assert(fov < 150);
		return radius / sin (fov / 360.0 * 3.14159) ;
	}
	
	void copy(const Camera& d)
	{
		_center         = d._center;
		_translate      = d._translate;
		_radius         = d._radius;
		_fov            = d._fov;
		_q              = d._q;
		return;
	}

};

#endif