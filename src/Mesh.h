/**
 * Mesh.h
 * (c)2008 Takashi Michikawa
 */
#ifndef __MESH_H__
#define __MESH_H__ 1
#include <OpenGL/gl.h>
#include <vecmath.h>
#include <deque>
#include <vector>
#include <math.h>
#include <iostream>

using std::cerr;
using std::operator<<;
using std::endl;

bool checkFormat(const char* filename, const char* format)
{
	return (strstr(filename, format) != NULL);
}

void getBoundingBox(std::deque<Vector3d>& vlist, Vector3d& bmin, Vector3d& bmax)
{
	bmin = vlist[0];
	bmax = vlist[0];
	for(size_t i = 1  ; i < vlist.size() ; i++)
	{
		if(vlist[i].x < bmin.x) bmin.x = vlist[i].x;
		if(vlist[i].y < bmin.y) bmin.y = vlist[i].y;
		if(vlist[i].z < bmin.z) bmin.z = vlist[i].z;
		if(vlist[i].x > bmax.x) bmax.x = vlist[i].x;
		if(vlist[i].y > bmax.y) bmax.y = vlist[i].y;
		if(vlist[i].z > bmax.z) bmax.z = vlist[i].z;
	}
	return;
}
/** 
 * @param[in] filename input obj file
 * @param[out] rad radis of bounding sphere
 * @param[out] center center of mass of obj file
 * @retval id object ID used in display list
 */
GLuint read_pointset(char* filename, double& rad, Vector3d& center )
{
	std::deque<Vector3d> vlist;
	
	cerr<<"Reading Point sets ";
	FILE *fp;
	if((fp = fopen(filename, "r")) == NULL )
	{
		cerr << "file cannot be read.\n";  exit(-1);
	}


	GLuint listid = glGenLists(1);
	glNewList(listid, GL_COMPILE);

	GLfloat mat_ambient[]    = { 0.2, 0.2, 0.2, 1.0 };
	GLfloat mat_diffuse[]    = { 0.8, 0.8, 0.8, 1.0 };
	GLfloat mat_specular[]   = { 0.8, 0.8, 1.0, 1.0 };
	GLfloat mat_shininess[]  = { 100.0 };

	glMaterialfv(GL_FRONT, GL_AMBIENT,  mat_ambient); 
	glMaterialfv(GL_FRONT, GL_DIFFUSE,  mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess);
	
	int num;
	fscanf(fp, "%d", &num);
	glBegin(GL_POINTS);
	glPointSize(0.1);
	for(int i = 0 ; i < num ; i++)
	{
		Vector3d p, n;
		fscanf(fp, "%lf %lf %lf", &(p.x), &(p.y), &(p.z));
		fscanf(fp, "%lf %lf %lf", &(n.x), &(n.y), &(n.z));
		vlist.push_back(p);
		n.normalize();
		glNormal3d(n.x,n.y,n.z);
		glVertex3d(p.x, p.y, p.z);
	}
	glEnd();
	fclose(fp);
	glEndList();

	
	Vector3d bmin,bmax;
	getBoundingBox(vlist, bmin, bmax);

	center =  0.5 * (bmin + bmax);

	Vector3d vec = bmax - bmin;
	rad = vec.length() * 0.5;
	cerr << "done...\n" << endl;
	return listid;
}
/** 
 * @param[in] filename input obj file
 * @param[out] rad radis of bounding sphere
 * @param[out] center center of mass of obj file
 * @retval id object ID used in display list
 */
GLuint read_mesh(char* filename, double& rad, Vector3d& center )
{
	std::deque<Vector3d> vlist;
	
	cerr<<"Reading Faces ";
	FILE *fp;
	char buf[512], dummy[32];

	if((fp = fopen(filename, "r")) == NULL )
	{
		cerr << "file cannot be read.\n";  exit(-1);
	}


	GLuint listid = glGenLists(1);
	glNewList(listid, GL_COMPILE);

	GLfloat mat_ambient[]    = { 0.2, 0.2, 0.2, 1.0 };
	GLfloat mat_diffuse[]    = { 0.8, 0.8, 0.8, 1.0 };
	GLfloat mat_specular[]   = { 0.8, 0.8, 1.0, 1.0 };
	GLfloat mat_shininess[]  = { 100.0 };

	glMaterialfv(GL_FRONT, GL_AMBIENT,  mat_ambient); 
	glMaterialfv(GL_FRONT, GL_DIFFUSE,  mat_diffuse);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS,mat_shininess);
	
	while( fgets(buf, 512, fp) )
	{
		if(  strncmp(buf, "#", 1) == 0 ) continue; //#‚Í–³Ž‹
		else if( strncmp(buf, "v", 1) == 0 )
		{	
			Vector3d p;
			sscanf(buf, "%s %lf %lf %lf", &dummy, &(p.x), &(p.y), &(p.z));
		
			vlist.push_back(p);
		}
		else if( strncmp(buf, "f", 1) ==0 )
		{
			std::vector<int> idx(3,0);
			sscanf(buf, "%s %d %d %d", &dummy, &(idx[0]), &(idx[1]), &(idx[2]));
			
	        glBegin(GL_POLYGON);
		    Vector3d p0 = vlist[idx[0] - 1];
			Vector3d p1 = vlist[idx[1] - 1];
			Vector3d p2 = vlist[idx[2] - 1];
			Vector3d v0 = p1- p0;
			Vector3d v1 = p2- p0;
			
			Vector3d n;
			n.cross(v0,v1);
			n.normalize();
			glNormal3d(n.x,n.y,n.z);
			glVertex3d(p0.x, p0.y, p0.z);
			glVertex3d(p1.x, p1.y, p1.z);
			glVertex3d(p2.x, p2.y, p2.z);
	        glEnd();
		}
	}
	fclose(fp);
    glEndList();

	
	Vector3d bmin,bmax;
	getBoundingBox(vlist, bmin, bmax);
	center =  0.5 * (bmin + bmax);

	Vector3d vec = bmax - bmin;
	rad = vec.length() * 0.5;
	cerr << "done...\n" << endl;
	return listid;
}
#endif
