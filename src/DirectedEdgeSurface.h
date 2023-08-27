///////////////////////////////////////////////////
//
//  Hamish Carr
//  September, 2020
//
//  ------------------------
//  DirectedEdgeSurface.h
//  ------------------------
//  
//  Base code for rendering assignments.
//
//  Minimalist (non-optimised) code for reading and 
//  rendering an object file
//  
//  We will make some hard assumptions about input file
//  quality. We will not check for manifoldness or 
//  normal direction, &c.  And if it doesn't work on 
//  all object files, that's fine.
//
//  While I could set it up to use QImage for textures,
//  I want this code to be reusable without Qt, so I 
//  shall make a hard assumption that textures are in 
//  ASCII PPM and use my own code to read them
//  
///////////////////////////////////////////////////

// include guard for DirectedEdgeSurface
#ifndef _DIRECTED_EDGE_SURFACE_H
#define _DIRECTED_EDGE_SURFACE_H

// include the C++ standard libraries we need for the header
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

#ifdef __APPLE__
#include <OpenGL/gl.h>
#else
#include <Windows.h>
#undef min
#undef max
#include <gl/GL.h>
#endif

// include the unit with Cartesian 3-vectors
#include "Cartesian3.h"
// the render parameters
#include "RenderParameters.h"
// the image class for a texture
#include "RGBAImage.h" 

// Half Edge Vertex Definition
struct HEVertex
{
    vec3 position;          // vertices positions
    unsigned int fd_edgeID; // first directed edgeIDs;

    HEVertex(vec3 _position, unsigned int _fd_edgeID)
    {
        position = _position;
        fd_edgeID = _fd_edgeID;
    }
    HEVertex(vec3 _position)
    {
        position = _position;
        fd_edgeID = -1;
    }
    HEVertex()
    {
        position = vec3();
        fd_edgeID = -1;
    }

};

// Half Edge Edge Definition
struct HEEdge
{
    unsigned int vert_ID;       // "Pointto" Vert ID;
    unsigned int face_ID;       // Adjcanit Face ID of this Edge;
    unsigned int next_edgeID;   // Next Half Edge ID of this edge;
    unsigned int prev_edgeID;   // Prev Edge ID of this edge;
    unsigned int oppo_edgeID;   // Opposite Edge ID of this edge;

    bool hasJoint;              // flag of this edge takes a joint vertex
    unsigned int jointID;       // actually is new vert ID

    HEEdge(unsigned int _vert, unsigned int _face, unsigned int _next, unsigned int _prev, unsigned int _oppo)
    {
        vert_ID = _vert;  face_ID = _face;

        next_edgeID = _next;  prev_edgeID = _prev; oppo_edgeID = _oppo;

        hasJoint = false; jointID = -1;

    }
    HEEdge()
    {
        vert_ID = -1; face_ID = -1;

        next_edgeID = -1; prev_edgeID = -1; oppo_edgeID = -1;

        hasJoint = false; jointID = -1;
    }
};

// Half Edge Mesh Definition
struct HalfEdgeMesh
{
    std::vector<HEVertex> verts;            // vertices data
    std::vector<HEEdge> edges;              // half-edges data
    std::vector<unsigned int> faces;        // faces data

    std::vector<vec3> normals;

    vec3 virutal_center;                   // center of the gravity
    float virutal_radius;                  // bounding sphere of the model

    //  constructor of mesh data structure
    HalfEdgeMesh()
    {
        verts.resize(0);
        edges.resize(0);

        normals.resize(0);
        faces.resize(0);

        virutal_center = vec3();
        virutal_radius = .0f;
    }
};

class DirectedEdgeSurface
{ // class DirectedEdgeSurface
public:

    const char* fileName = nullptr;

    HalfEdgeMesh mesh;

    // constructor will initialise to safe values
    DirectedEdgeSurface();
    
    // process of loop subdivision
    void loop_subdivision();

    // read routine returns true on success, failure otherwise
    bool ReadObjectStream(std::istream &geometryStream);

    // write routine
    void WriteObjectStream(std::ostream &geometryStream);

    // routine to render
    void Render(RenderParameters *renderParameters);

    void Write2System()
    {
        if(fileName!=nullptr)
        {
            std::string fileString(fileName);
            std::string output_filename = fileString.substr(fileString.find_last_of("/\\") + 1);

            std::ofstream outfile("./" + output_filename , std::ofstream::binary);

            WriteObjectStream(outfile);
        }
    }

private:
    // get neighbors around a vertex. It should return its neighbors and the degree of this vertex
    std::vector<unsigned int> umbrella(HEVertex vert, unsigned int vertID, HalfEdgeMesh mesh);

}; // class DirectedEdgeSurface

// end of include guard for DirectedEdgeSurface
#endif
