#pragma once

#include "utils.h"
#include "triangle.h"

std::list<Triangle> find_triangles(Graph &g);

std::list<Triangle> find_triangles_naive(const Graph &g);

