// SPDX-License-Identifier: GPL-3.0-or-later
// CLI adapter for the pinned fork; local changes are documented beside CMakeLists.txt.
#include <AutoRemesher/AutoRemesher>
#define TINYOBJLOADER_IMPLEMENTATION
#include <tiny_obj_loader.h>
#include <chrono>
#ifdef _WIN32
#include <codecvt>
#include <locale>
#endif
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>

struct Progress {
    bool quiet = false;
    std::mutex mutex;
    int last = -1;
    std::string status;
};

static void progress(void* tag, float fraction, const char* message)
{
    auto& state = *static_cast<Progress*>(tag);
    if (state.quiet) return;
    std::lock_guard<std::mutex> lock(state.mutex);
    const int percent = std::max(0, std::min(100, int(fraction * 100)));
    const std::string status = message ? message : "Remeshing";
    if (percent == state.last && status == state.status) return;
    state.last = percent;
    state.status = status;
    std::cout << "PROGRESS " << percent << ' ' << status << std::endl;
}

static double number(const std::string& text, double low, double high)
{
    size_t consumed = 0;
    double value = std::stod(text, &consumed);
    if (consumed != text.size() || !std::isfinite(value) || value < low || value > high)
        throw std::runtime_error("Option out of range: " + text);
    return value;
}

// Keep Windows paths in their native UTF-16 representation. UTF-8 conversion is
// only for option parsing and diagnostics; filenames never take an ANSI detour.
#ifdef _WIN32
using NativePath = std::wstring;
static std::string argumentText(const wchar_t* text)
{
    return std::wstring_convert<std::codecvt_utf8_utf16<wchar_t>>().to_bytes(text);
}
int wmain(int argc, wchar_t** argv)
#else
using NativePath = std::string;
static std::string argumentText(const char* text) { return text; }
int main(int argc, char** argv)
#endif
{
    try {
        NativePath input, output, report;
        std::string model = "organic";
        size_t quads = 5000;
        bool clothFolds = true;
        double scaling = 0, adaptivity = 1, anisotropy = 1, sharp = 90, smooth = 0;
        Progress state;
        for (int i = 1; i < argc; ++i) {
            std::string option = argumentText(argv[i]);
            if (option == "--version") {
                std::cout << "autoremesher-fine-details " << ENGINE_REVISION << " rebuilt-cli-2 " << ENGINE_PATCH << "\n";
                return 0;
            }
            if (option == "--help") {
                std::cout << "AutoRemesherEngine --input in.obj --output out.obj\n"
                    "  --quads N --scaling F --adaptivity F --anisotropy F\n"
                    "  --sharp-edge-degrees D --smooth-normal-degrees D\n"
                    "  --cloth-folds (default) | --no-cloth-folds\n"
                    "  --model-type organic|hard --report out.json --quiet --version\n";
                return 0;
            }
            if (option == "--quiet") { state.quiet = true; continue; }
            if (option == "--cloth-folds") { clothFolds = true; continue; }
            if (option == "--no-cloth-folds") { clothFolds = false; continue; }
            if (++i >= argc) throw std::runtime_error("Missing value for " + option);
            if (option == "--input") { input = argv[i]; continue; }
            if (option == "--output") { output = argv[i]; continue; }
            if (option == "--report") { report = argv[i]; continue; }
            std::string value = argumentText(argv[i]);
            if (option == "--model-type") model = value;
            else if (option == "--quads") {
                double n = number(value, 1, 100000000);
                if (std::floor(n) != n) throw std::runtime_error("Quad count must be an integer");
                quads = size_t(n);
            }
            else if (option == "--scaling") scaling = number(value, .25, 4);
            else if (option == "--adaptivity") adaptivity = number(value, 0, 1);
            else if (option == "--anisotropy") anisotropy = number(value, 0, 1);
            else if (option == "--sharp-edge-degrees") sharp = number(value, 0, 180);
            else if (option == "--smooth-normal-degrees") smooth = number(value, 0, 180);
            else throw std::runtime_error("Unknown option: " + option);
        }
        if (input.empty() || output.empty() || input == output)
            throw std::runtime_error("Distinct --input and --output paths are required");
        if (model != "organic" && model != "hard") throw std::runtime_error("Invalid model type");
        tinyobj::attrib_t attributes;
        std::vector<tinyobj::shape_t> shapes;
        std::vector<tinyobj::material_t> materials;
        std::string warning, error;
        // Windows C++ runtimes provide native wchar_t filename overloads in
        // C++14. The parser and its defaults stay the same as the filename API.
        std::ifstream inputFile(input.c_str());
        if (!inputFile)
            throw std::runtime_error("OBJ load failed: Cannot open file [" + argumentText(input.c_str()) + "]");
        // The prior LoadObj call supplied no material base directory: preserve
        // MaterialFileReader's working-directory lookup and filename semantics.
        tinyobj::MaterialFileReader materialReader("");
        if (!tinyobj::LoadObj(&attributes, &shapes, &materials, &warning, &error,
                &inputFile, &materialReader))
            throw std::runtime_error("OBJ load failed: " + error);
        inputFile.close();
        std::vector<AutoRemesher::Vector3> vertices;
        std::vector<std::vector<size_t>> triangles;
        for (size_t i = 0; i < attributes.vertices.size(); i += 3) {
            for (size_t axis = 0; axis < 3; ++axis)
                if (!std::isfinite(attributes.vertices[i + axis])) throw std::runtime_error("Non-finite input vertex");
            vertices.emplace_back(attributes.vertices[i], attributes.vertices[i + 1], attributes.vertices[i + 2]);
        }
        for (const auto& shape : shapes) {
            if (shape.mesh.indices.size() % 3) throw std::runtime_error("OBJ triangulation failed");
            for (size_t i = 0; i < shape.mesh.indices.size(); i += 3) {
                std::vector<size_t> face;
                for (size_t k = 0; k < 3; ++k) {
                    const int v = shape.mesh.indices[i + k].vertex_index;
                    if (v < 0 || size_t(v) >= vertices.size()) throw std::runtime_error("Invalid input face index");
                    face.push_back(size_t(v));
                }
                triangles.push_back(face);
            }
        }
        if (triangles.empty()) throw std::runtime_error("Empty input mesh");
        const auto start = std::chrono::steady_clock::now();
        AutoRemesher::AutoRemesher engine(vertices, triangles);
        engine.setTargetTriangleCount(quads * 2);
        if (scaling > 0) engine.setScaling(scaling);
        engine.setGradientAdaptivity(adaptivity);
        engine.setAnisotropy(anisotropy);
        engine.setClothFoldGuidance(clothFolds);
        engine.setSharpEdgeDegrees(sharp);
        engine.setSmoothNormalDegrees(smooth);
        engine.setModelType(model == "hard" ? AutoRemesher::ModelType::HardSurface : AutoRemesher::ModelType::Organic);
        engine.setTag(&state);
        engine.setProgressHandler(progress);
        if (!engine.remesh() || engine.remeshedQuads().empty()) throw std::runtime_error("Remesher returned no mesh");
        size_t quadCount = 0;
        for (const auto& p : engine.remeshedVertices())
            for (size_t axis = 0; axis < 3; ++axis)
                if (!std::isfinite(p[axis])) throw std::runtime_error("Non-finite output vertex");
        for (const auto& face : engine.remeshedQuads()) {
            if (face.size() < 3) throw std::runtime_error("Invalid output polygon");
            for (size_t v : face)
                if (v >= engine.remeshedVertices().size()) throw std::runtime_error("Invalid output face index");
            quadCount += face.size() == 4;
        }
        std::ofstream out(output.c_str());
        out << std::setprecision(17);
        for (const auto& p : engine.remeshedVertices()) out << "v " << p.x() << ' ' << p.y() << ' ' << p.z() << '\n';
        for (const auto& face : engine.remeshedQuads()) {
            out << 'f';
            for (size_t v : face) out << ' ' << v + 1;
            out << '\n';
        }
        out.close();
        if (!out) throw std::runtime_error("OBJ write failed");
        const double seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
        const std::string summary = "{\"engine_revision\":\"" ENGINE_REVISION "\",\"quads\":"
            + std::to_string(quadCount) + ",\"faces\":" + std::to_string(engine.remeshedQuads().size())
            + ",\"cloth_folds\":" + (clothFolds ? "true" : "false")
            + ",\"vertices\":" + std::to_string(engine.remeshedVertices().size())
            + ",\"seconds\":" + std::to_string(seconds) + "}";
        std::cout << summary << std::endl;
        for (const auto& line : engine.phaseReport()) std::cerr << line << '\n';
        if (!report.empty()) {
            std::ofstream file(report.c_str());
            file << summary << '\n';
            file.close();
            if (!file) throw std::runtime_error("Report write failed");
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "AutoRemesherEngine: " << error.what() << std::endl;
        return 2;
    }
}
