#include "DelauneyTriangulationDivConq.h"

#include <random>

int main(int argc, char** argv)
{
    std::vector<double2> points;

    const auto n = 20;
    const double A = 10.;

    std::random_device rd;
    std::mt19937 gen(0);

    std::uniform_real_distribution<> dis(-A, A);

    points.reserve(n);
    for (size_t i = 0; i < n; ++i)
    {
        double2 const pos{ dis(gen), dis(gen)};
        points.emplace_back(pos);
    }

	DelauneyTriangulationDivConq::TriangulationDelauneyBuilder builder;
    try
    {
        auto triangulation = builder.build(points);

        triangulation.print("graph");
    }
    catch (const std::runtime_error& e)
    {
        std::cerr << "Could not build, error: " << e.what() << std::endl;
    }

	return 0;
}