#ifndef BOUNDINGBOX_H
#define BOUNDINGBOX_H
#include <array>

class BoundingBox {

    private:
        std::array<float, 3> min;
        std::array<float, 3> max;

    public:

        BoundingBox() : min({0, 0, 0}), max({0, 0, 0}) {}

        BoundingBox(const std::array<float, 3> &min, const std::array<float, 3> &max) : min(min), max(max) {}

        bool intersects(const Ray &ray) const {
            float tmin = std::numeric_limits<float>::lowest();
            float tmax = std::numeric_limits<float>::max();

            for (int i = 0; i < 3; i++) {
                float t1 = (min[i] - ray.origin()[i]) / ray.direction()[i];
                float t2 = (max[i] - ray.origin()[i]) / ray.direction()[i];

                if (t1 > t2) std::swap(t1, t2);

                tmin = std::max(tmin, t1);
                tmax = std::min(tmax, t2);

                if (tmax < tmin) return false;
            }

            return true;
        }

        void expand(const BoundingBox &other) {
            for (int i = 0; i < 3; i++) {
                min[i] = std::min(min[i], other.min[i]);
                max[i] = std::max(max[i], other.max[i]);
            }
        }

        void draw(const BoundingBox &box) {
            std::array<float, 24> vertices = {
                box.min[0], box.min[1], box.min[2],
                box.max[0], box.min[1], box.min[2],
                box.max[0], box.max[1], box.min[2],
                box.min[0], box.max[1], box.min[2],
                box.min[0], box.min[1], box.max[2],
                box.max[0], box.min[1], box.max[2],
                box.max[0], box.max[1], box.max[2],
                box.min[0], box.max[1], box.max[2]
            };

            std::array<unsigned int, 36> indices = {
                0, 1, 2, 2, 3, 0,
                4, 5, 6, 6, 7, 4,
                4, 5, 1, 1, 0, 4,
                6, 7, 3, 3, 2, 6,
                1, 5, 6, 6, 2, 1,
                4, 0, 3, 3, 7, 4
            };

            glBegin(GL_LINE_STRIP);
            for (unsigned int i : indices) {
                glVertex3d(vertices[3 * i], vertices[3 * i + 1], vertices[3 * i + 2]);
            }
            glEnd();
        }


        std::pair<std::array<float, 3>, std::array<float, 3>> getBounds(const Mesh &mesh) {
            std::array<float, 3> min = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
            std::array<float, 3> max = {std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest(), std::numeric_limits<float>::lowest()};

            for (unsigned int i : mesh.triangles_array) {
                const MeshVertex &vertex = mesh.vertices[i];
                for (int j = 0; j < 3; j++) {
                    min[j] = std::min(min[j], vertex.position[j]);
                    max[j] = std::max(max[j], vertex.position[j]);
                }
            }

            return {min, max};
        }
};

#endif