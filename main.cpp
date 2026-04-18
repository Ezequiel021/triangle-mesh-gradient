#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <map>
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
};

struct edge
{
    int a;
    int b;
};

int load_vtk()
{
    std::ifstream vtk("mesh.vtk");
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

    std::vector<float2> vertices(vertex_count);

    for (int i = 0; i < vertex_count; i++)
    {
        vtk >> vertices[i].x;
        vtk >> vertices[i].y;
        vtk >> discard;
    }

    // Skip "CELLS" line
    std::getline(vtk, discard);
    std::vector<polygon> faces;
    std::map<int64_t, int32_t> mp;

    while (! vtk.eof())
    {
        polygon poly_buffer;
        vtk >> poly_buffer.vertex_count;
        poly_buffer.vertex_id.resize(vertex_count);
        
        for (int i : poly_buffer.vertex_id)
        {
            vtk >> i;
        }
    }
}

int main()
{
    int32_t a = 20232321;
    int32_t b = 12312123;
    std::cout << std::bitset<64>(hash(a, b)) << "\n";
    std::cout << std::bitset<64>(hash(b, a)) << "\n";
    std::cout << std::bitset<32>(a) << std::bitset<32>(b) << "\n";
}