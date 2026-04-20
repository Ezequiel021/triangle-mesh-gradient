# triangle-mesh-gradient

https://github.com/Ezequiel021/triangle-mesh-gradient

# Dependencies

OpenMPI y CMake

# Build

```
git clone https://github.com/Ezequiel021/triangle-mesh-gradient && cd triangle-mesh-gradient
mkdir build && cd build
cmake build .
```

# Run
```
mpirun -np number-of-process ./gradiente filename.vtk
```

Reemplaze filename.vtk con el archivo vtk que quiera procesar. El repositorio incluye la malla de mantarraya y una malla de ejemplo:

# Example
```
mpirun -np 4 ./gradiente ../mallaraya.vtk
```