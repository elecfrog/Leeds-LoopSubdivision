///////////////////////////////////////////////////
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

// include the header file
#include "DirectedEdgeSurface.h"

// include the C++ standard libraries we want
#include <iostream>
#include <string>

// include the Cartesian 3- vector class
#include "Cartesian3.h"
#include "SphereVertices.h"

#define MAXIMUM_LINE_LENGTH 1024

std::vector<unsigned int> DirectedEdgeSurface::umbrella(HEVertex vert, unsigned int vertID, HalfEdgeMesh mesh)
{
    //unsigned int degree = 0;
    std::vector<unsigned int> neighbors;
    unsigned int origin_edgeID = vert.fd_edgeID;
    unsigned int origin_vertID = vertID;
    unsigned int des_edgeID = origin_edgeID;

//    printf("origin_edgeID %d\t, origin_vertID %d\n", origin_edgeID, origin_vertID);

    do
    {
//        printf("edgeID - %d | vert ID - %d | oppo_edgeID - %d | ", des_edgeID, mesh.edges[des_edgeID].vert_ID, mesh.edges[des_edgeID].oppo_edgeID);

        neighbors.push_back(mesh.edges[des_edgeID].vert_ID);

        //degree++;

        des_edgeID = mesh.edges[des_edgeID].oppo_edgeID; //get the oppo edge ID, using oppo back home

        if (mesh.edges[des_edgeID].vert_ID == origin_vertID) // check it if back home
        {
//            printf("next_edgeID %d\n", mesh.edges[des_edgeID].next_edgeID);

            des_edgeID = mesh.edges[des_edgeID].next_edgeID; // get the next edge ID
        }
        else { break; }


    } while (des_edgeID != origin_edgeID);

    //printf("degree %d\n", degree);
    //mesh.verts[vertID].degree = degree;
//    for (auto& v : neighbors)
//        printf("vertID %d\t", v);
//    printf("\n");

    //mesh.verts[vertID].neighbors = vert_list;
    return neighbors;
}

// constructor will initialise to safe values
DirectedEdgeSurface::DirectedEdgeSurface()
{ // DirectedEdgeSurface()

} // DirectedEdgeSurface()

// read routine returns true on success, failure otherwise
bool DirectedEdgeSurface::ReadObjectStream(std::istream &geometryStream)
    { // ReadObjectStream()
    
    // create a read buffer
    char readBuffer[MAXIMUM_LINE_LENGTH];
    
    // the rest of this is a loop reading lines & adding them in appropriate places
    while (true)
        { // not eof
		// token for identifying meaning of line
		std::string token;

        // character to read
        geometryStream >> token;
        
        // check for eof() in case we've run out
        if (geometryStream.eof())
            break;

        // otherwise, switch on the token we read
		if (token == "#")
			{ // comment 
			// read and discard the line
			geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
            } // comment
		else if (token == "Vertex")
			{ // vertex
			// variables for the read
			unsigned int vertexID;
			geometryStream >> vertexID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
            if (vertexID != mesh.verts.size())
				{ // bad vertex ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad vertex ID				
			
			// read in the new vertex position
			Cartesian3 newVertex;
			geometryStream >> newVertex;
            HEVertex new_vert = HEVertex(newVertex);
			// and add it to the vertices
            mesh.verts.push_back(new_vert);
			} // vertex
		else if (token == "Normal")
			{ // normal
			// variables for the read
			unsigned int normalID;
			geometryStream >> normalID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
            if (normalID != mesh.normals.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
			// read in the new normal
			Cartesian3 newNormal;
			geometryStream >> newNormal;
			
			// and add it to the vertices
            mesh.normals.push_back(newNormal);
			} // normal
		else if (token == "FirstDirectedEdge")
			{ // first directed edge
			// variables for the read
			unsigned int FDEID;
			geometryStream >> FDEID;
//            std::cout<<FDEID<<std::endl;

            // it has to be next valid 0-based ID, so
			// reject line if it isn't
            if (mesh.verts[FDEID].fd_edgeID != -1)
            { // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
            } // bad ID
			
			// read in the new FDE
			unsigned int newFDE;
			geometryStream >> newFDE;
//            std::cout<<newFDE<<std::endl;

            // and add it to the vertices
            mesh.verts[FDEID].fd_edgeID = newFDE;

        } // first directed edge
		else if (token == "Face")
			{ // face
			// variables for the read
			unsigned int faceID;
			geometryStream >> faceID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
            if (faceID != mesh.faces.size()/3)
				{ // bad face ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad face ID				
			
			// read in the new face vertex (3 times)
			unsigned int newFaceVertex;
			geometryStream >> newFaceVertex;
            mesh.faces.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
            mesh.faces.push_back(newFaceVertex);
			geometryStream >> newFaceVertex;
            mesh.faces.push_back(newFaceVertex);
			} // face
		else if (token == "OtherHalf")
			{ // other half
			// variables for the read
			unsigned int otherHalfID;
			geometryStream >> otherHalfID;
			// it has to be next valid 0-based ID, so
			// reject line if it isn't
            if (otherHalfID != mesh.edges.size())
				{ // bad ID
				// read and discard the line
				geometryStream.getline(readBuffer, MAXIMUM_LINE_LENGTH);
				} // bad ID				
			
            // Do my own version here:
			unsigned int newOtherHalf;
			geometryStream >> newOtherHalf;


            unsigned int faceID = otherHalfID/3;
            unsigned int next;
            unsigned int prev;
            unsigned int vertID = mesh.faces[3*faceID + otherHalfID%3];
            if (otherHalfID % 3 == 0) // current edge is e0 ,next is e1, prev is e2
            {
                next = otherHalfID + 1;
                prev = otherHalfID + 2;
            }
            else if (otherHalfID % 3 == 1) // current edge is e1, next is e2, prev is e0
            {
                next = otherHalfID + 1;
                prev = otherHalfID - 1;
            }
            else if (otherHalfID % 3 == 2) // current edge is e2, next is e0, prev is e1
            {
                next = otherHalfID - 2; // e0
                prev = otherHalfID - 1; // e1
            }

            // Vert ID | Face ID | Next ID | Prev ID | Oppo ID
            HEEdge new_edge = HEEdge(vertID, faceID, next, prev, newOtherHalf);
            mesh.edges.push_back(new_edge);
			} // other half
        } // not eof

    // if there are any vertices at all
    if (mesh.verts.size() != 0)
    { // non-empty vertex set
        // sum up all of the vertex positions

        for (unsigned int vertex = 0; vertex < mesh.verts.size(); vertex++)
            mesh.virutal_center = mesh.virutal_center + mesh.verts[vertex].position;

        // and divide through by the number to get the average position
        // also known as the barycentre
        mesh.virutal_center = mesh.virutal_center / mesh.verts.size();

        if(int(mesh.virutal_center.x) == int(mesh.virutal_center.y) == int(mesh.virutal_center.z) == 0)
            mesh.virutal_center = vec3();

        // now compute the largest distance from the origin to a vertex
        for (unsigned int vertex = 0; vertex < mesh.verts.size(); vertex++)
            { // per vertex
            // compute the distance from the barycentre
            float distance = (mesh.verts[vertex].position - mesh.virutal_center).length();
            
            // now test for maximality
            if (distance > mesh.virutal_radius)
                mesh.virutal_radius = distance;
            } // per vertex
        } // non-empty vertex set

    // return a success code

    /* PRINT INFO OF READED DATA!*/
//    for(int v = 0; v< mesh.verts.size(); ++v)
//    {
//        std::cout<< v << " | "<< mesh.verts[v].position << std::endl;
//    }
//    for(int v = 0; v< mesh.verts.size(); ++v)
//    {
//        std::cout<< v << " | "<< mesh.verts[v].fd_edgeID<< std::endl;
//    }
//    for(int fde = 0; fde< mesh.fd_edges.size(); ++fde)
//        std::cout<< fde << "  "<< mesh.fd_edges[fde] <<std::endl;
//    for(int n = 0; n< mesh.normals.size(); ++n)
//        std::cout<< n << "  "<< mesh.normals[n] <<std::endl;
//    for(int f = 0; f< mesh.faces.size()/3; ++f)
//        std::cout<< f << "  "<< mesh.faces[3*f] <<" "<< mesh.faces[3*f+1]<<" "<< mesh.faces[3*f+2]<<std::endl;
//    for(int oh = 0; oh< mesh.otherhalf.size(); ++oh)
//        std::cout<< oh << "  "<< mesh.otherhalf[oh] <<std::endl;
//    std::cout<< "virutal_center" << mesh.virutal_center <<std::endl;
//    std::cout<< "virutal_radius" << mesh.virutal_radius <<std::endl;

    return true;
    } // ReadObjectStream()

// write routine
void DirectedEdgeSurface::WriteObjectStream(std::ostream &geometryStream)
    { // WriteObjectStream()
	geometryStream << "#" << std::endl; 
    geometryStream << "# Loop Subdivision" << std::endl;
    geometryStream << "# Fanxiang Zhou" << std::endl;
    geometryStream << "# 201598922" << std::endl;
    geometryStream << "# Surface vertices=" << mesh.verts.size() << " faces=" << mesh.faces.size()/3 << std::endl;
	geometryStream << "#" << std::endl; 

	// output the vertices
    for (unsigned int vertex = 0; vertex < mesh.verts.size(); vertex++)
        geometryStream << "Vertex " << vertex << " " << std::fixed << mesh.verts[vertex].position << std::endl;

    // and the normal vectors
    for (unsigned int normal = 0; normal < mesh.normals.size(); normal++)
        geometryStream << "Normal " << normal << " " << std::fixed << mesh.normals[normal] << std::endl;

	// and the first directed edges
    for (unsigned int vertex = 0; vertex < mesh.verts.size(); vertex++)
        geometryStream << "FirstDirectedEdge " << vertex<< " " << std::fixed << mesh.verts[vertex].fd_edgeID << std::endl;

    // and the faces - increment is taken care of internally
    for (unsigned int face = 0; face < mesh.faces.size(); )
        { // per face
        geometryStream << "Face " << face << " ";
        
        // read in three vertices
        geometryStream << mesh.faces[face++] << " ";
        geometryStream << mesh.faces[face++] << " ";
        geometryStream << mesh.faces[face++];
            
        geometryStream << std::endl;
        } // per face

	// and the other halves
    for (unsigned int dirEdge = 0; dirEdge < mesh.edges.size(); dirEdge++)
        geometryStream << "OtherHalf " << dirEdge << " " << mesh.edges[dirEdge].oppo_edgeID << std::endl;
} // WriteObjectStream()

// routine to render
void DirectedEdgeSurface::Render(RenderParameters *renderParameters)
    { // Render()
    // Ideally, we would apply a global transformation to the object, but sadly that breaks down
    // when we want to scale things, as unless we normalise the normal vectors, we end up affecting
    // the illumination.  Known solutions include:
    // 1.   Normalising the normal vectors
    // 2.   Explicitly dividing the normal vectors by the scale to balance
    // 3.   Scaling only the vertex position (slower, but safer)
    // 4.   Not allowing spatial zoom (note: sniper scopes are a modified projection matrix)
    //
    // Inside a game engine, zoom usually doesn't apply. Normalisation of normal vectors is expensive,
    // so we will choose option 2.  

    // Scale defaults to the zoom setting
    float scale = renderParameters->zoomScale;
    scale /= mesh.virutal_radius;
        
    //  now scale everything
    glScalef(scale, scale, scale);

    // apply the translation to the centre of the object if requested
    glTranslatef(-mesh.virutal_center.x, -mesh.virutal_center.y, -mesh.virutal_center.z);

    // start rendering
    glBegin(GL_TRIANGLES);

	// set colour for pick render - ignored for regular render
	glColor3f(1.0, 1.0, 1.0);

    // loop through the faces
    for (unsigned int face = 0; face < mesh.faces.size(); face +=3)
		{ // per face
		// if we want flat normals, compute them here
		if (renderParameters->useFlatNormals)
			{ // flat normals
			// find two vectors along edges of the triangle
            Cartesian3 pq = mesh.verts[mesh.faces[face+1]].position - mesh.verts[mesh.faces[face]].position;
            Cartesian3 pr = mesh.verts[mesh.faces[face+2]].position - mesh.verts[mesh.faces[face]].position;

			// take their cross product and normalise
			Cartesian3 faceNormal = pq.cross(pr).unit();

			// and use it to set the glNormal
			glNormal3f(faceNormal.x * scale, faceNormal.y * scale, faceNormal.z * scale);
			} // flat normals

		// we have made a HARD assumption that we have enough normals
		for (unsigned int vertex = face; vertex < face+3; vertex++)
			{ // per vertex
		
			// if we are using smooth normals
			if (!renderParameters->useFlatNormals)
				// set the normal vector
				glNormal3f
					(
                    mesh.normals[mesh.faces[vertex]].x * scale,
                    mesh.normals[mesh.faces[vertex]].y * scale,
                    mesh.normals[mesh.faces[vertex]].z * scale
					);
			
			// and set the vertex position
			glVertex3f
				(
                mesh.verts[mesh.faces[vertex]].position.x,
                mesh.verts[mesh.faces[vertex]].position.y,
                mesh.verts[mesh.faces[vertex]].position.z
				);

			} // per vertex

		} // per face

    // close off the triangles
    glEnd();
    
    // now we add a second loop to render the vertices if desired
    if (!renderParameters->showVertices)
    	return;

	glDisable(GL_LIGHTING);

	// loop through the vertices
    for (unsigned int vertex = 0; vertex < mesh.verts.size(); vertex++)
		{ // per vertex
		// use modelview matrix (not most efficient solution, but quickest to code)
		glMatrixMode(GL_MODELVIEW);
		glPushMatrix();
        glTranslatef(mesh.verts[vertex].position.x, mesh.verts[vertex].position.y, mesh.verts[vertex].position.z);
		glScalef(0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize, 0.1 * renderParameters->vertexSize);
		renderTriangulatedSphere();
		glPopMatrix();
		} // per vertex 
    
} // Render()

void DirectedEdgeSurface::loop_subdivision()
{
    // 1. Create New Joints (vetices)
    unsigned int joint_counts = 0;
    for (auto& edge : mesh.edges)
    {
        if (edge.hasJoint == false)
        {
            unsigned int curr_vertID = edge.vert_ID;
            unsigned int oppo_vertID = this->mesh.edges[edge.oppo_edgeID].vert_ID;
            unsigned int next_vertID = this->mesh.edges[edge.next_edgeID].vert_ID;
            unsigned int opxt_vertID = this->mesh.edges[this->mesh.edges[edge.oppo_edgeID].next_edgeID].vert_ID;
            //  edge & oppo_edge handle this new vert
            vec3 joint_position =
                (mesh.verts[curr_vertID].position + mesh.verts[oppo_vertID].position) * 0.375f +
                (mesh.verts[next_vertID].position + mesh.verts[opxt_vertID].position) * 0.125f ;

            // insert normals at the same time
            vec3 new_normal = vec3();

            mesh.normals.push_back(new_normal);

            HEVertex new_joint = HEVertex(joint_position); // new a joint vert and insert into the vert list

            edge.jointID = mesh.verts.size();
            mesh.edges[edge.oppo_edgeID].jointID = edge.jointID;

            edge.hasJoint = true;
            mesh.edges[edge.oppo_edgeID].hasJoint = true;

            mesh.verts.push_back(new_joint);
            joint_counts++;
        }
    }

    unsigned int ID = 0;
    for (auto& i : mesh.edges)
        printf("e%d| f%d| prev %d| v%d| next %d| oppo %d | jointID %d \n",
            ID++, i.face_ID, i.prev_edgeID, i.vert_ID, i.next_edgeID, i.oppo_edgeID, i.jointID);

    // 2. Update Old Verts : Each Vert use umbrella
    ID = 0;
    for (auto& vert : mesh.verts)
    {
        if (ID < mesh.verts.size() - joint_counts)
        {
            std::vector<unsigned int> neigbors = umbrella(vert, ID++, this->mesh);
            // n = degree
            unsigned int degree = neigbors.size();
            float u = 3.0f / float(8 * degree);
            //printf("u %f", u);
            if (degree == 3) u = 0.1875;

            vec3 sum_of_neigbors = vec3();
            for (auto& neigbor : neigbors)
            {
                sum_of_neigbors = sum_of_neigbors + mesh.verts[neigbor].position;
            }

            vert.position = (1.0f - float(degree) * u) * vert.position + u * sum_of_neigbors;
        }

    }

    // 3. Redefine each face f_k will be subdivided in 4
    unsigned int current_edge_size = mesh.edges.size();
    for (unsigned int f = 0; f < current_edge_size / 3; ++f)
    {
        // get 6 vert IDs from origin faces and edges
        unsigned int v_x = mesh.faces[3 * f];                // 0
        unsigned int v_y = mesh.faces[3 * f + 1];            // 1
        unsigned int v_z = mesh.faces[3 * f + 2];            // 2

        unsigned int v_u = mesh.edges[3 * f].jointID;        // 12
        unsigned int v_v = mesh.edges[3 * f + 1] .jointID;   // 13
        unsigned int v_w = mesh.edges[3 * f + 2 ].jointID;   // 14

        /*
                 v_y
                /   \
          e_2 v_w -- v_v  e_1
              /  \  /  \
           v_z -- v_u -- v_x
                  e_0
        */

        // the central face retaining the original index
        // original face
        mesh.faces[3 * f] = v_u;     // 12
        mesh.faces[3 * f + 1] = v_v; // 13
        mesh.faces[3 * f + 2] = v_w; // 14

        // new face f_aplha  0-13-12
        mesh.faces.push_back(v_x); mesh.faces.push_back(v_v);mesh.faces.push_back(v_u);

        // new face f_theta 1-14-13
        mesh.faces.push_back(v_y); mesh.faces.push_back(v_w);mesh.faces.push_back(v_v);

        // new face f_gamma 2-12-14
        mesh.faces.push_back(v_z); mesh.faces.push_back(v_u);mesh.faces.push_back(v_w);
    }
    for (unsigned int f = 0; f < mesh.faces.size() / 3; ++f)
        printf("faces %d | %d %d %d \n", f, mesh.faces[3 * f], mesh.faces[3 * f + 1], mesh.faces[3 * f + 2]);

    // 4. Rebuild half edges
    mesh.edges.resize(0);
    for (unsigned int f = 0; f < mesh.faces.size() / 3; ++f)
    {
        //EDGE      VERT_ID                 FACEID              NEXT_ID     PREV_ID                      OPPO_ID(with default 0)
        HEEdge e0(mesh.faces[3 * f],        f, 3 * f + 1, 3 * f + 2,  0);
        HEEdge e1(mesh.faces[3 * f + 1],    f, 3 * f + 2, 3 * f    ,  0);
        HEEdge e2(mesh.faces[3 * f + 2],    f, 3 * f    , 3 * f + 1,  0);
        mesh.edges.push_back(e0);
        mesh.edges.push_back(e1);
        mesh.edges.push_back(e2);
    }

    for (unsigned int index = 0; index < mesh.edges.size() ; ++index)
    {
        // v_start
        unsigned int v_start = mesh.edges[mesh.edges[index].prev_edgeID].vert_ID;

        // v_point
        unsigned int v_point = mesh.edges[index].vert_ID;

        //printf("start point %d, point to %d\n", v_start, v_point);

        for (unsigned int inner = 0; inner < mesh.edges.size(); ++inner)
        {
            unsigned int i_start = mesh.edges[mesh.edges[inner].prev_edgeID].vert_ID;
            if (i_start == v_point && mesh.edges[inner].vert_ID == v_start)
            {
                mesh.edges[index].oppo_edgeID = inner;
                break;

            }
        }
    }

    ID = 0;
    for (auto& edge : mesh.edges)
    printf("e%d| v%d| f%d| next %d| prev %d| oppo %d \n",
        ID++,edge.vert_ID, edge.face_ID, edge.next_edgeID, edge.prev_edgeID, edge.oppo_edgeID);

    // 5. update normals update_normal();

    // 6. Update FirstDirectedEdge
    for (unsigned int vertID = 0;  vertID < mesh.verts.size(); ++vertID)
    {
        bool stop = false;
        for (unsigned int edgeID = 0; edgeID < mesh.edges.size(); ++edgeID)
        {
            unsigned int start_vert = mesh.edges[mesh.edges[edgeID].prev_edgeID].vert_ID;
            if (start_vert == vertID && stop == false)
            {
                mesh.verts[vertID].fd_edgeID = edgeID;
                stop = true;
            }
        }
    }

    ID = 0;
    for(int v = 0; v< mesh.verts.size(); ++v)
    {
        std::cout<< v << " | "<< mesh.verts[v].position << std::endl;
    }
}
