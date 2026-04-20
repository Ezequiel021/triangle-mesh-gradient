#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <bitset>
#include <mpi.h>
#include <string>
#include <cmath>
#define EPS 10e-15
// bit interleaving unordered pair hashing
inline int64_t hash(int32_t a, int32_t b)
{
    return a < b ? (int64_t)a ^ ((int64_t)b << 32) : (int64_t)b ^ ((int64_t)a << 32);
}

struct float2
{
    double y;
    double x;

    float2() {}
    float2(double x, double y)
    {
        this->x = x;
        this->y = y;
    }

    void set(double x, double y)
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
    double f; ///Valor f en el centroide;
    float2 c;///Centroide
    //Derivadas que calcularemos de centroide.
    double grad_x;
    double grad_y;
};

struct edge
{
    int a;
    int b;
};
///Funcion auxiliar calcula el centroide dado una cara
float2 centroid(const polygon& p, const std::vector<float2>& vertices)
{
    double x = 0.0, y = 0.0;
    for (int v : p.vertex_id){
        x += vertices[v].x;
        y += vertices[v].y;
    }
    return float2(x / p.vertex_count, y / p.vertex_count);
}

double funct(double x,double y){
    return 20.0*(std::sin(4*x)*std::sin(3*y)+0.3*std::cos(6*x)*std::sin(5*y)+1);
}

int load_vtk(const std::string namefile,std::vector<float2> &vertices, std::vector<polygon> &faces)
{
    std::ifstream vtk(namefile);
    if (!vtk.is_open())
    {
        return -1;
    }

    std::string buffer;
    ///Skipeamos header
    std::getline(vtk,buffer);
    std::getline(vtk,buffer);
    std::getline(vtk,buffer);
    std::getline(vtk,buffer);
    int vertex_c;
    //input es de la forma  points <vertex_c> <type>
    vtk >> buffer >> vertex_c >> buffer;

    vertices.resize(vertex_c);

    ///Leemos vertices
    for(int i = 0; i < vertex_c;i++){
        double x,y,z;
        vtk >> x >> y >> z;
        vertices[i].set(x,y);
    }

    ///Leemos las celdas (caras)
    int num_cells,total_size;
    vtk >> buffer >> num_cells >> total_size;
    faces.clear();
    faces.reserve(num_cells);

    for(int i = 0; i < num_cells;i++){
        polygon poly;
        vtk >> poly.vertex_count;
        poly.vertex_id.resize(poly.vertex_count);
        for(int j = 0; j < poly.vertex_count;j++){
            vtk >> poly.vertex_id[j];
        }
        if(poly.vertex_count == 3)//Ignoramos las aristas, solo triangulos
            faces.push_back(poly);
    }
    std::cout << "Loaded VTK: " << vertex_c << " vertices, " << num_cells << " cells\n";
    return 0;
}

void find_neighbor_poly(std::vector<polygon> &faces)
{
    std::unordered_map<int64_t, int32_t> mp;

    int face_count = faces.size();

    for (int f = 0; f < face_count; f++)
    {
        int n = faces[f].vertex_count;

        for (int i = 0; i < n; i++)
        {
            int a = faces[f].vertex_id[i];
            int b = faces[f].vertex_id[(i + 1) % n];

            int64_t key = hash(a, b);

            auto it = mp.find(key);

            if (it != mp.end())
            {
                int other = it->second;

                faces[f].neighboor_id.push_back(other);
                faces[other].neighboor_id.push_back(f);

                mp.erase(it);
            }
            else
            {
                mp[key] = f;
            }
        }
    }
}
//Funcion auxiliar, le da a cada cara su centroide, y f evaluado en ese centroide.
void setFacesVal(std::vector<polygon>& faces, std::vector<float2>&vertices){
    for(int i = 0; i < faces.size();i++){
        faces[i].c = centroid(faces[i],vertices);
        double x = faces[i].c.x;
        double y = faces[i].c.y;
        faces[i].f = funct(x,y);
    }
}

//Dado un triangulo, calcula el gradiente en el centroide
void compute_grad(std::vector<polygon>& faces){
    int N  = faces.size();

    for(int i = 0; i < N;i++){
        double x0 = faces[i].c.x;
        double y0 = faces[i].c.y;
        double f0 = faces[i].f;
        double Sxx = 0.0; double Syy = 0.0;
        double Sxy = 0.0;
        double Sxf = 0.0; double Syf = 0.0;
        int M = faces[i].neighboor_id.size();
        for(int k = 0; k < M;k++){
            int j = faces[i].neighboor_id[k];
            double dx = faces[j].c.x - x0;
            double dy = faces[j].c.y - y0;
            double df = faces[j].f - f0;
            Sxx += dx * dx;
            Syy += dy * dy;
            Sxy += dx * dy;

            Sxf += dx * df;
            Syf += dy * df;
        }
        double D = Sxx*Syy-Sxy*Sxy;
        faces[i].grad_x = (Syy * Sxf - Sxy * Syf) / D;
        faces[i].grad_y = (Sxx * Syf - Sxy * Sxf) / D;

    }
}
void save_vertices(
    const std::vector<float2>& vertices,
    const std::vector<double>& node_grad_x,
    const std::vector<double>& node_grad_y,
    const std::vector<double>& exact_grad_x,
    const std::vector<double>& exact_grad_y,
    const std::vector<double>& error_x,
    const std::vector<double>& error_y)
{
    std::ofstream file("output.vtk");

    file << "# vtk DataFile Version 2.0\n";
    file << "Gradiente en vtx\n";
    file << "ASCII\n";
    file << "DATASET POLYDATA\n";

    file << "POINTS " << vertices.size() << " double\n";
    for (auto& v : vertices)
        file << v.x << " " << v.y << " 0.0\n";

    file << "VERTICES " << vertices.size() << " " << 2 * vertices.size() << "\n";
    for (int i = 0; i < (int)vertices.size(); i++)
        file << "1 " << i << "\n";

    file << "POINT_DATA " << vertices.size() << "\n";

    // Guardamos gradiente aproximado
    file << "VECTORS grad_approx double\n";
    for (int i = 0; i < (int)vertices.size(); i++)
        file << node_grad_x[i] << " " << node_grad_y[i] << " 0.0\n";

    // Guardamos gradiente analitico
    file << "VECTORS grad_exact double\n";
    for (int i = 0; i < (int)vertices.size(); i++)
        file << exact_grad_x[i] << " " << exact_grad_y[i] << " 0.0\n";

    // Error absoluto
    file << "VECTORS error_vec double\n";
    for (int i = 0; i < (int)vertices.size(); i++)
        file << error_x[i] << " " << error_y[i] << " 0.0\n";

    file.close();
}

int main()
{
    std::vector<float2> vertices;
    std::vector<polygon> faces;

    load_vtk("Malla_Mantarraya.vtk",vertices, faces);

    find_neighbor_poly(faces);
    setFacesVal(faces,vertices);
    compute_grad(faces);

    ///Para el paso de interpolacion, primero hacemos vector, tal que\
    dado el idx de un vertice, tendremos idx de las caras que lo tienen.
    std::vector<std::vector<int>> vert_to_face(vertices.size());
    for(int i = 0; i < faces.size(); i++){
        for(int v : faces[i].vertex_id){
            vert_to_face[v].push_back(i);
        }
    }
    std::vector<double> node_grad_x(vertices.size(), 0.0);
    std::vector<double> node_grad_y(vertices.size(), 0.0);

    ///El paso de interpolacion, promediamos gradientes de los centroides a triangulos\
    vecinos del vertice
    for(int v = 0; v < vertices.size(); v++){
        for(int fi : vert_to_face[v]){
            node_grad_x[v] += faces[fi].grad_x;
            node_grad_y[v] += faces[fi].grad_y;
        }
        node_grad_x[v] /= vert_to_face[v].size();
        node_grad_y[v] /= vert_to_face[v].size();
    }

    std::vector<double> exact_grad_x(vertices.size());
    std::vector<double> exact_grad_y(vertices.size());
    std::vector<double> error_x(vertices.size());
    std::vector<double> error_y(vertices.size());

    ///Calculamos error abs de solucion analitica vs aproximada
    for(int v = 0; v < (int)vertices.size(); v++){
        double x = vertices[v].x;
        double y = vertices[v].y;
        exact_grad_x[v] = 20.0*(4.0*std::cos(4*x)*std::sin(3*y) - 1.8*std::sin(6*x)*std::sin(5*y));
        exact_grad_y[v] = 20.0*(3.0*std::sin(4*x)*std::cos(3*y) + 1.5*std::cos(6*x)*std::cos(5*y));
        error_x[v] = std::abs(node_grad_x[v] - exact_grad_x[v]);
        error_y[v] = std::abs(node_grad_y[v] - exact_grad_y[v]);
    }

    ///Guardamos para cada vertice, los gradientes x,y analiticos, y los aproximados, con el error.
    save_vertices(vertices, node_grad_x, node_grad_y, exact_grad_x, exact_grad_y, error_x, error_y);
}
