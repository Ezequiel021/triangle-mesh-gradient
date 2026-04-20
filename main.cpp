#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include <bitset>
#include <mpi.h>
#include <chrono>
#include <cmath>
#include <algorithm>

// Intercalación de bits para hashing de pares no ordenados
inline int64_t hash(int32_t a, int32_t b)
{
    return a < b ? (int64_t)a ^ ((int64_t)b << 32) : (int64_t)b ^ ((int64_t)a << 32);
}

struct vec2
{
    double x; // o double y, según tu orden
    double y;

    // Inicializar obligatoriamente en 0
    vec2() : x(0.0), y(0.0) {}
    vec2(double x, double y) : x(x), y(y) {}
};

struct polygon
{
    int vertex_count;
    std::vector<int> vertex_id;
    std::vector<int> neighboor_id;
};

// Función analítica f(x,y)
double f_analytical(double x, double y)
{
    return 20.0 * (std::sin(4.0 * x) * std::sin(3.0 * y) + 0.3 * std::cos(6.0 * x) * std::sin(5.0 * y) + 1.0);
}

// Derivadas parciales analíticas para el cálculo del error
vec2 grad_analytical(double x, double y)
{
    double df_dx = 20.0 * (4.0 * std::cos(4.0 * x) * std::sin(3.0 * y) - 1.8 * std::sin(6.0 * x) * std::sin(5.0 * y));
    double df_dy = 20.0 * (3.0 * std::sin(4.0 * x) * std::cos(3.0 * y) + 1.5 * std::cos(6.0 * x) * std::cos(5.0 * y));
    return vec2(df_dx, df_dy);
}

    // ...
int load_vtk(const std::string &filename, std::vector<vec2> &vertices, std::vector<polygon> &faces)
{
    std::ifstream vtk(filename);
    if (!vtk.is_open())
        return -1;

    std::string buffer, discard;
    while (vtk >> buffer && buffer != "POINTS")
        ;

    int vertex_count;
    vtk >> vertex_count >> discard;

    vertices.resize(vertex_count);
    for (int i = 0; i < vertex_count; i++)
    {
        float z; // Descartar Z ya que es 2D
        vtk >> vertices[i].x >> vertices[i].y >> z;
    }

    while (vtk >> buffer && buffer != "CELLS")
    {
    }

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

        if (poly_buffer.vertex_count == 3)
        {
            faces.push_back(poly_buffer);
        }
    }
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
                faces[f].neighboor_id.push_back(nf->second);
                faces[nf->second].neighboor_id.push_back(f);
                mp.erase(nf);
            }
            else
            {
                mp.insert({key, f});
            }
        }
    }
}

// Genera el archivo VTK con los resultados nodales
void export_vtk(const std::string &filename, const std::vector<vec2> &vertices, const std::vector<polygon> &faces,
                const std::vector<vec2> &approx_grad, const std::vector<vec2> &exact_grad)
{
    std::ofstream out(filename);
    out << "# vtk DataFile Version 3.0\nGradient Results\nASCII\nDATASET UNSTRUCTURED_GRID\n";

    out << "POINTS " << vertices.size() << " double\n";
    for (const auto &v : vertices)
    {
        out << v.x << " " << v.y << " 0.0\n";
    }

    int total_indices = 0;
    for (const auto &f : faces)
        total_indices += f.vertex_count + 1;

    out << "\nCELLS " << faces.size() << " " << total_indices << "\n";
    for (const auto &f : faces)
    {
        out << f.vertex_count;
        for (int vid : f.vertex_id)
            out << " " << vid;
        out << "\n";
    }

    out << "\nCELL_TYPES " << faces.size() << "\n";
    for (size_t i = 0; i < faces.size(); ++i)
        out << "5\n"; // 5 = VTK_TRIANGLE

    out << "\nPOINT_DATA " << vertices.size() << "\n";

    // Campo: Gradiente Aproximado
    out << "VECTORS Gradiente_Aproximado double\n";
    for (const auto &g : approx_grad)
        out << g.x << " " << g.y << " 0.0\n";

    // Campo: Gradiente Analitico
    out << "VECTORS Gradiente_Analitico double\n";
    for (const auto &g : exact_grad)
        out << g.x << " " << g.y << " 0.0\n";

    // Campo: Error Absoluto
    out << "VECTORS Error_Absoluto double\n";
    for (size_t i = 0; i < vertices.size(); ++i)
    {
        out << std::abs(approx_grad[i].x - exact_grad[i].x) << " "
            << std::abs(approx_grad[i].y - exact_grad[i].y) << " 0.0\n";
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<vec2> vertices;
    std::vector<polygon> faces;

    // Todos leen la malla para tener acceso implícito a "ghost nodes"
    if (load_vtk(argv[1], vertices, faces) != 0)
    {
        if (rank == 0)
            std::cerr << "Error cargando la malla.\n";
        MPI_Finalize();
        return -1;
    }

    find_neighbor_poly(faces);

    int total_faces = faces.size();
    std::vector<vec2> centroids(total_faces);
    std::vector<double> f_vals(total_faces);

    // Calcular centroides y valores de la función
    for (int i = 0; i < total_faces; ++i)
    {
        double cx = 0, cy = 0;
        for (int vid : faces[i].vertex_id)
        {
            cx += vertices[vid].x;
            cy += vertices[vid].y;
        }
        cx /= faces[i].vertex_count;
        cy /= faces[i].vertex_count;
        centroids[i] = vec2(cx, cy);
        f_vals[i] = f_analytical(cx, cy);
    }

    // Partición del dominio (elementos)
    int local_n = total_faces / size;
    int remainder = total_faces % size;
    int start_idx = rank * local_n + std::min(rank, remainder);
    int end_idx = start_idx + local_n + (rank < remainder ? 1 : 0);
    int local_count = end_idx - start_idx;

    std::vector<vec2> local_gradients(local_count);

    MPI_Barrier(MPI_COMM_WORLD);
    auto start_time = std::chrono::high_resolution_clock::now();

    // Cálculo del gradiente por mínimos cuadrados en paralelo
    for (int i = start_idx; i < end_idx; ++i)
    {
        // Justo antes de calcular los mínimos cuadrados para cada i:
        if (faces[i].vertex_count != 3)
        {
            local_gradients[i - start_idx] = vec2(0.0f, 0.0f);
            continue;
        }
        double x0 = centroids[i].x;
        double y0 = centroids[i].y;
        double f0 = f_vals[i];

        double Sxx = 0, Syy = 0, Sxy = 0, Sxf = 0, Syf = 0;

        for (int neighbor_idx : faces[i].neighboor_id)
        {
            double dx = centroids[neighbor_idx].x - x0;
            double dy = centroids[neighbor_idx].y - y0;
            double df = f_vals[neighbor_idx] - f0;

            Sxx += dx * dx;
            Syy += dy * dy;
            Sxy += dx * dy;
            Sxf += dx * df;
            Syf += dy * df;
        }

        double D = Sxx * Syy - Sxy * Sxy;
        double a = 0, b = 0;

        // Evitar división por cero en casos degenerados
        if (std::abs(D) > 1e-15)
        {
            a = (Syy * Sxf - Sxy * Syf) / D;
            b = (Sxx * Syf - Sxy * Sxf) / D;
        }

        local_gradients[i - start_idx] = vec2(a, b);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    auto end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end_time - start_time;

    // Comunicación Punto a Punto: Enviar resultados al maestro
    std::vector<vec2> global_gradients;

    if (rank == 0)
    {
        global_gradients.resize(total_faces);
        // Copiar datos locales del maestro
        for (int i = 0; i < local_count; ++i)
        {
            global_gradients[start_idx + i] = local_gradients[i];
        }

        // Recibir datos del resto de procesos (Punto a Punto)
        for (int p = 1; p < size; ++p)
        {
            int p_start = p * local_n + std::min(p, remainder);
            int p_count = local_n + (p < remainder ? 1 : 0);
            if (p_count > 0)
            {
                MPI_Recv(&global_gradients[p_start], p_count * 2, MPI_DOUBLE, p, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }
    }
    else
    {
        if (local_count > 0)
        {
            MPI_Send(local_gradients.data(), local_count * 2, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
        }
    }

    // Interpolación a nodos y escritura (Solo Maestro)
    if (rank == 0)
    {
        std::cout << "Tiempo de calculo paralelo (" << size << " hilos): " << diff.count() << " s\n";

        std::vector<vec2> node_grad_approx(vertices.size(), vec2(0, 0));
        std::vector<int> node_element_count(vertices.size(), 0);

        // Interpolar centroides a nodos
        for (int i = 0; i < total_faces; ++i)
        {
            for (int vid : faces[i].vertex_id)
            {
                node_grad_approx[vid].x += global_gradients[i].x;
                node_grad_approx[vid].y += global_gradients[i].y;
                node_element_count[vid]++;
            }
        }

        std::vector<vec2> node_grad_exact(vertices.size());
        for (size_t v = 0; v < vertices.size(); ++v)
        {
            if (node_element_count[v] > 0)
            {
                node_grad_approx[v].x /= node_element_count[v];
                node_grad_approx[v].y /= node_element_count[v];
            }
            node_grad_exact[v] = grad_analytical(vertices[v].x, vertices[v].y);
        }

        export_vtk("resultados.vtk", vertices, faces, node_grad_approx, node_grad_exact);
        std::cout << "Resultados guardados en resultados.vtk listo para ParaView.\n";
    }

    MPI_Finalize();
    return 0;
}