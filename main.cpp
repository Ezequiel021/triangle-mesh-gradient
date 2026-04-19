#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <bitset>
#include <mpi.h>
#include <chrono>

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
    std::ifstream vtk("../mallaraya.vtk");
    if (!vtk.is_open())
    {
        std::cerr << "Error: No se pudo abrir ../mesh.vtk\n";
        return -1;
    }

    std::string buffer, discard;

    // Skip headers til keyword "POINTS"
    while (vtk >> buffer && buffer != "POINTS");

    int vertex_count;
    vtk >> vertex_count;
    vtk >> discard;

    vertices.resize(vertex_count);
    for (int i = 0; i < vertex_count; i++)
    {
        vtk >> vertices[i].x;
        vtk >> vertices[i].y;
        vtk >> discard;
    }

    // Skip til kw "CELLS"
    while (vtk >> buffer && buffer != "CELLS"){}

    int cell_count, cell_list_size;
    vtk >> cell_count >> cell_list_size;

    faces.reserve(cell_count);
    for (int c = 0; c < cell_count; c++)
    {
        polygon poly_buffer;
        vtk >> poly_buffer.vertex_count;
        poly_buffer.vertex_id.resize(poly_buffer.vertex_count);

        for (int i = 0; i < poly_buffer.vertex_count; i++)
        {
            vtk >> poly_buffer.vertex_id[i];
        }

        faces.push_back(poly_buffer);
    }

    std::cout << "Archivo .vtk cargado con exito. Vertices: "
              << vertex_count << ", Caras: " << cell_count << "\n";
    return 0;
}
void find_neighbor_poly(std::vector<polygon> &faces)
{
    std::unordered_map<int64_t, int32_t> mp;
    int face_count = faces.size();
    for (int f = 0; f < face_count; f++)
    {
        int v_count = faces[f].vertex_id.size();
        for (int i = 0; i < v_count; i++)
        {
            int32_t v1 = faces[f].vertex_id[i];
            int32_t v2 = faces[f].vertex_id[(i + 1) % v_count];

            int64_t key = hash(v1, v2);
            
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

int main()
{
    std::vector<float2> vertices;
    std::vector<polygon> faces;

    load_vtk(vertices, faces);

    auto start = std::chrono::high_resolution_clock::now();
    find_neighbor_poly(faces);
    auto end = std::chrono::high_resolution_clock::now();

    auto runtime_neighbors = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    std::cout << "Tiempo en calcular los polígonos vecinos: " << runtime_neighbors.count() << "\n";
}