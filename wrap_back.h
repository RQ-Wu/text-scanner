#pragma once
#include"global.h"

Mesh init_mesh(int rows, int cols);
void wrap_mesh(Mesh& mesh, DisMap displacementMap, Mat mask);
void fixMesh(Mat mask, meshPos& meshVertexCoord);