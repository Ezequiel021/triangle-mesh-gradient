#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <bitset>
#include <mpi.h>

// bit interleaving unordered pair hashing
inline int64_t hash(int32_t a, int32_t b)
{
    return a < b ? (int64_t)a ^ ((int64_t)b << 32) : (int64_t)b ^ ((int64_t)a << 32);
}

struct float2
{
    float y;
    float x;

    float2() {}
    float2(float x, float y)
    {
        this->x = x;
        this->y = y;
    }

    void set(float x, float y)
    {    
        this->x = x;
        this->y = y;
    }
};

struct polygon
{
    int vertex_count;
    std::vector<int> vertex_id;
    std::vector<int> neighboor_id;
};

struct edge
{
    int a;
    int b;
};

int load_vtk(std::vector<float2> &vertices, std::vector<polygon> &faces)
{
    std::ifstream vtk("../mesh.vtk");
    if (!vtk.is_open())
    {
        return -1;
    }

    std::string buffer, discard;

    // Skip "ASCII" line
    std::getline(vtk, discard);
    
    // Skip "DATASET ... " line
    std::getline(vtk, discard);

    int vertex_count;
    // Read "POINTS ... " line
    vtk >> discard;
    vtk >> vertex_count;
    vtk >> discard;

    vertices.resize(vertex_count);
    for (int i = 0; i < vertex_count; i++)
    {
        vtk >> vertices[i].x;
        vtk >> vertices[i].y;
        vtk >> discard;
    }

    // Skip "CELLS" line
    std::getline(vtk, discard);

    std::string sbuf;
    while (std::getline(vtk, sbuf))
    {
        std::stringstream ss(sbuf);

        polygon poly_buffer;
        ss >> poly_buffer.vertex_count;
        poly_buffer.vertex_id.resize(vertex_count);
        
        for (int i : poly_buffer.vertex_id)
        {
            ss >> i;
        }

        faces.push_back(poly_buffer);
    }

    std::cout << "archivo .vtk cargado con exito" << "\n";
    return 0;
}

void find_neighbor_poly(std::vector<polygon> &faces)
{
    std::unordered_map<int64_t, int32_t> mp;
    int face_count = faces.size();
    for (int f = 0; f < face_count; f++)
    {
        for (int i = 0; i < faces[f].vertex_count; i++)
        {
            for (int j = i + 1; j < faces[f].vertex_count; j++)
            {
                int64_t key = hash(faces[f].vertex_id[i], faces[f].vertex_id[j]);
                auto nf = mp.find(key);
                if (nf != mp.end())
                {
                    faces[f].neighboor_id.emplace_back(nf->second);
                    faces[nf->second].neighboor_id.push_back(f); // Ensure symmetry
                    mp.erase(nf);
                }
                else
                {
                    mp.insert({key, f});
                }
            }
        }
    }
}

int main()
{
    std::vector<float2> vertices;
    std::vector<polygon> faces;

    load_vtk(vertices, faces);

    find_neighbor_poly(faces);
    for (int i = 0; i < faces.size(); i++)
    {
        std::cout << i << ": ";
        for (int j = 0; j < faces[i].neighboor_id.size(); j++)
        {
            std::cout << faces[i].neighboor_id[j] << ", ";
        }
        std::cout << "\n";
    }
}