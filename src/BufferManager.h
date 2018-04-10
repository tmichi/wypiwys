/**
* BufferManager.h
* (c)2008 Takashi Michikawa
*/
#ifndef __BUFFER_MANAGER_H__
#define __BUFFER_MANAGER_H__ 1

#include <OpenGL/gl.h>
#include <vector>
#include <vecmath.h>
using namespace kh_vecmath;


class BufferManager
{
private:
	int _width;
	int _height;
	std::vector<GLfloat> _depth; //depth buffer 
	std::vector<bool> _mask;    //maskimage of painting
	GLdouble _modelview[16];
	GLdouble _projection[16];
public:

	BufferManager()
	{ 

	};
	~BufferManager(){};


	bool startPaintMode()
	{
		glDisable(GL_LIGHTING);

		int viewport[4];
		glGetIntegerv(GL_VIEWPORT,viewport);

		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluOrtho2D(0,viewport[2],0,viewport[3]);

		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		return true;
	}

	bool startViewMode()
	{
		glEnable(GL_LIGHTING);	
		glDepthFunc(GL_LESS);
		return true;
	}

	bool readDepth() 
	{
		_depth.clear();
		int viewport[4];
		glGetIntegerv(GL_VIEWPORT,viewport); 
		_width  = viewport[2];
		_height = viewport[3];

		glGetDoublev(GL_MODELVIEW_MATRIX, _modelview);
		glGetDoublev(GL_PROJECTION_MATRIX, _projection);

		const int npixel = _width * _height;
		_depth.assign(npixel, 0.0f);
		glReadBuffer(GL_FRONT);
		glReadPixels(0, 0, _width, _height, GL_DEPTH_COMPONENT, GL_FLOAT, &(_depth[0]));
		glClear(GL_DEPTH_BUFFER_BIT);

		return true;
	}
	bool readMask() 
	{
		int npixel = _width * _height;
		_mask.clear();
		_mask.assign(npixel, false);

		std::vector<GLfloat> buf(npixel, 0.0f);
		glReadBuffer(GL_FRONT);
		glReadPixels(0, 0, _width, _height, GL_DEPTH_COMPONENT, GL_FLOAT, &(buf[0])); 

		for(size_t i = 0 ; i < buf.size() ; i++)
		{	
			if(0 < buf[i] && buf[i] < 1)_mask[i] = true;
		}
		return true;
	}

	bool getPoints(std::vector<Vector3d>& pnt, std::vector<double>& weight)
	{
		pnt.clear();
		weight.clear();


		int viewport[4];
		glGetIntegerv(GL_VIEWPORT,viewport); 


		std::vector<Vector3d> pointArray;
		//points
		for(size_t i = 0 ; i < _depth.size() ; i++)
		{	
			if(_depth[i] <= 0 || 1 <= _depth[i] ) _mask[i] = false;		
		}

		for(int h = 0 ; h < _height ; h++)
		{
			for(int w = 0 ; w < _width ; w++)
			{	       
				Vector3d p;
				int ww = w;
				int hh =  h ;
				int id =  h * _width + w;

				gluUnProject((double)ww, (double)hh , _depth[id],
					_modelview, _projection, viewport,
					&(p.x), &(p.y), &(p.z) );
				pointArray.push_back(p);
			}
		}


		for(int h = 0 ; h < _height - 1 ; h++)
		{
			for(int w = 0 ; w < _width - 1 ; w++)
			{	       
				std::vector<int> nbrid;
				int dh[] = {0, 1, 1, 0};
				int dw[] = {0, 0, 1, 1};

				int count = 0;
				for(size_t i = 0 ; i < 4 ; i++)
				{					
					int ww = w + dw[i];
					int hh = h + dh[i];			
					nbrid.push_back(ww + hh  * _width);

					if(_mask[nbrid[i] ])
					{
						count++;
					}
				}

				if(count > 3)
				{
					Vector3d center;
					for(size_t i = 0 ; i < 4 ; i++)
					{
						if(_mask[nbrid[i] ] )center += pointArray[nbrid[i]];
					}
					center.scale( 1.0 / count);

					nbrid.push_back(nbrid[0]); // to avoid using %
					double area = 0 ; 
					for(size_t i= 0 ; i < 4 ; i++)
					{
						if(_mask[nbrid[i]] && _mask[nbrid[i+1]])
						{
							Vector3d v0 = pointArray[nbrid[i]  ] - center;
							Vector3d v1 = pointArray[nbrid[i+1]] - center;
							Vector3d n;
							n.cross(v0,v1);
							area += n.length();
						}
					}
					pnt.push_back(center);
					weight.push_back(area);
				}
			}
		}
		pointArray.clear();
		_mask.clear();
		_depth.clear();
		return true;
	};
};

bool dumpImage()
{
	static int count = 0;
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT,viewport); 
	int width  = viewport[2];
	int height = viewport[3];
	const int npixel = width * height;

	std::vector<GLubyte> rgb(npixel*3, 0x00);
	glReadBuffer(GL_FRONT);
	glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, &(rgb[0]));

	char file[256];
	sprintf(file, "dump-%d.raw", count);
	count++;
	FILE* fp = fopen(file,"wb"); 
	fwrite(&(rgb[0]), 1, 3*npixel, fp);
	fclose(fp);
	return true;
}

bool dumpDepth()
{
	static int count2 = 0;
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT,viewport); 
	int width  = viewport[2];
	int height = viewport[3];
	const int npixel = width * height;

	std::vector<GLfloat> buf(npixel, 0.0f);
	glReadBuffer(GL_FRONT);
	glReadPixels(0, 0, width, height, GL_DEPTH_COMPONENT, GL_FLOAT, &(buf[0]));

	char file[256];
	sprintf(file, "depth-%d.raw", count2);
	count2++;
	FILE* fp = fopen(file, "wb"); 
	for(size_t t = 0 ; t < buf.size() ; t++)
	{	
		unsigned short v = (unsigned short)(buf[t] * 65535) ; 
		fwrite(&v, 2, 1, fp);
	}
	fclose(fp);
	return true;
}


#endif // __BUFFER_MANAGER_H__
