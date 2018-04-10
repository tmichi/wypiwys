/**
 * CameraAnimation.h
 * (c)2008 Takashi Michikawa
 */
#ifndef __CAMERA_ANIMATION_H__
#define __CAMERA_ANIMATION_H__ 1
#include "Camera.h"
#include <vector>
class CameraAnimation
{
private:
	double   _dRadius;
	double   _dFOV;
	Vector3d _dCenter;
	Quat4d   _dQ;

	Camera* _camera;

	int _max_step;
	int _step;
public:
	CameraAnimation()
	{
		_max_step = -1;	
		_step = 0;
	};
	~CameraAnimation(){};

	void init(Camera* curCamera, const std::vector<Vector3d> &pnt, const std::vector<double> &weight, const int maxstep)
	{
		_camera = curCamera;
		if(maxstep <= 0 ) _max_step = 100;	
		_max_step = maxstep;
		_step = 0;

		const Camera newCamera = this->getNewCamera(_camera, pnt, weight);
		
		double subdiv = 1.0 / _max_step;
		_dCenter = (newCamera.getCenter() - curCamera->getCenter() ) * subdiv;
		_dRadius = (newCamera.getRadius() - curCamera->getRadius() ) * subdiv;
		_dFOV    = (newCamera.getFOV()    - curCamera->getFOV()    ) * subdiv;
		_dQ = newCamera.getQuaternion();
		_dQ.mulInverse(curCamera->getQuaternion());
		_dQ.interpolate(Quat4d(0,0,0,1), 1.0 - subdiv);
	};

	bool animate()
	{
		if(_step <= _max_step )
		{
			_step++;
			_camera->move(_dRadius, _dCenter, _dQ, _dFOV);
			_camera->updatePerspective();
			return true;
		}
		return false;
	}

private:
	// implementation of wypiwys
	Camera getNewCamera(const Camera* curCamera, const std::vector<Vector3d> &pnt, const std::vector<double> &weight ) const
	{
		const Vector3d curRay = curCamera->getRay();

		Vector3d newCenter;
		double sumWeight = 0.0;
		for(size_t i = 0 ; i < weight.size() ; i++)
		{
			newCenter += weight.at(i) * pnt.at(i);
			sumWeight += weight.at(i) ;
		}
		newCenter.scale(1.0 / sumWeight) ;

		//compute radius
		double newRadius = 0.0 ;
		for(size_t i = 0 ; i < weight.size() ; i++)
		{
			const Vector3d v = pnt.at(i) - newCenter;
			if(newRadius < v.length()) newRadius = v.length();
		}

		//compute rotation quaternion		
		const Vector3d newRay = this->computeRay(curRay, newCenter, pnt, weight);
		double angle = curRay.angle(newRay);
		Quat4d newq(0.0,0.0,0.0,1.0);
		if(angle > 1.0e-10)
		{
			Vector3d axis, right;
			axis.cross(curRay, newRay);
			right.cross(axis,  curRay);
			if(right.dot(newRay) < 0 ) axis.scale(-1);//angle *= -1;
			newq.set(AxisAngle4d(axis, angle));	
		}
		newq.mul(curCamera->getQuaternion());

		return Camera(newRadius, newCenter, newq, curCamera->getFOV());

	}
	Vector3d 
		get_eigenvector(const double e, const Matrix3d& m) const
	{	
		Vector3d eigenvector;
		const double a = 1;
		const double b = -m.m00 - m.m11 - m.m22;
		const double c = -m.m01 * m.m10 - m.m20 * m.m02 - m.m12 * m.m21
			+m.m00 * m.m11 + m.m11 * m.m22 + m.m00 * m.m22;
/*		const double d = 
			-m.m00 * m.m11 * m.m22 + m.m00 * m.m12 * m.m21 
			-m.m01 * m.m12 * m.m20 + m.m01 * m.m10 * m.m22
			-m.m02 * m.m10 * m.m21 + m.m02 * m.m11 * m.m20;
*/
		const double bb = b + a * e;
		const double cc = c + bb * e; 

		Matrix3d e0;
		e0.setIdentity();
		Matrix3d a2 = m * m + bb * m + cc * e0;

		eigenvector.set(a2.m00, a2.m10, a2.m20);
		if(eigenvector.lengthSquared() < 1.0e-16)
		{
			eigenvector.set(a2.m01, a2.m11, a2.m21);
			if(eigenvector.lengthSquared() < 1.0e-16)
			{
				eigenvector.set(a2.m02, a2.m12, a2.m22);
				if(eigenvector.lengthSquared() < 1.0e-16)return eigenvector; 
			}
		}
		eigenvector.normalize();
		return eigenvector;    
	}

	std::vector<double> 
		get_eigenvalue(const Matrix3d& m) const
	{
		std::vector<double> eigenvalue;

		const double a = 1;
		const double b = -m.m00 - m.m11 - m.m22;
		const double c = -m.m01 * m.m10 - m.m20 * m.m02 - m.m12 * m.m21
			+m.m00 * m.m11 + m.m11 * m.m22 + m.m00 * m.m22;
		const double d = 
			-m.m00 * m.m11 * m.m22 + m.m00 * m.m12 * m.m21 
			-m.m01 * m.m12 * m.m20 + m.m01 * m.m10 * m.m22
			-m.m02 * m.m10 * m.m21 + m.m02 * m.m11 * m.m20;

		double x1 = 0;
		double x2 = 1;

		if (d != 0)
		{     
			while(fabs(x1-x2) > 1.0e-10)
			{
				x2=x1-(((a*x1*x1*x1)+(b*x1*x1)+(c*x1)+d)/((3*a*x1*x1)+(2*b*x1)+c));
				x1=x2;
			}
		}
		x2 = 0;

		double aa = a;   
		double bb = b + a * x1;
		double cc = c + bb * x1; 

		eigenvalue.push_back(x1);
		double D = bb * bb - 4 * aa * cc;

		if(fabs(D) < 1.0e-15 )
		{
			eigenvalue.push_back( -bb / ( 2 * aa));
		}
		else if(D>0)
		{
			eigenvalue.push_back( (-bb - sqrt(D)) / (2 * aa) );
			eigenvalue.push_back( (-bb + sqrt(D)) / (2 * aa) );
		}
		std::sort(eigenvalue.begin(), eigenvalue.end());
		return eigenvalue;
	}
	Vector3d 
		computeRay(const Vector3d& curRay, const Vector3d& center, const std::vector<Vector3d>& pnt, const std::vector<double>& weight)const
	{
		// covariance matrix
		Matrix3d m;
		for(size_t i = 0 ; i < pnt.size() ; i++)
		{
			const Vector3d v = pnt.at(i) - center;
			Matrix3d dm(v.x * v.x , v.x * v.y, v.x * v.z, 
				v.y * v.x , v.y * v.y, v.y * v.z, 
				v.z * v.x , v.z * v.y, v.z * v.z);
			dm.mul(weight.at(i));
			m.add(dm);
		}

		const std::vector<double> eigenvalue = get_eigenvalue(m);
		Vector3d newRay = get_eigenvector(eigenvalue.at(0), m);
		if( newRay.dot(curRay) < 0 ) newRay.scale(-1);
		newRay.normalize();
		return newRay;
	}
};
#endif //__CAMERA_ANIMATION_H__