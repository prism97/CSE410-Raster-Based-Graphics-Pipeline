#include <bits/stdc++.h>
#include "bitmap_image.hpp"

using namespace std;

#define PI 3.14159265358979323846

double eyeX, eyeY, eyeZ;
double lookX, lookY, lookZ;
double upX, upY, upZ;
double fovY, aspectRatio, near, far;

int screenWidth, screenHeight;
double xLimit, yLimit;
double zFrontLimit, zRearLimit;

ofstream stageFile;

class Color {
public:
    int r, g, b;
    Color();
    Color(int r, int g, int b);
};

Color::Color() {

}

Color::Color(int r, int g, int b) {
    this->r = r;
    this->g = g;
    this->b = b;
}

class Point {
public:
    double x, y, z, w;
    Point();
    Point(double x, double y, double z);
    Point(const Point &p1);
    double getCoordinate(int index) const;
    void setCoordinate(int index, double val);
    void wScale();
    void print() const;
    double dotProduct(const Point& p) const;
    Point crossProduct(const Point& p) const;
    Point scalarMultiply(double s) const;
    Point subtract(const Point& p) const;
    void normalize();
};

Point::Point() {
    this->w = 1;
}

Point::Point(double x, double y, double z) {
    this->x = x;
    this->y = y;
    this->z = z;
    this->w = 1;
}

Point::Point(const Point &p1) {
    this->x = p1.x;
    this->y = p1.y;
    this->z = p1.z;
    this->w = p1.w;
}

double Point::getCoordinate(int index) const {
    switch (index) {
        case 0: return this->x;
        case 1: return this->y;
        case 2: return this->z;
        case 3: return this->w;
        default: return 0.0;
    }
}

void Point::setCoordinate(int index, double val) {
    switch (index) {
        case 0: this->x = val; break;
        case 1: this->y = val; break;
        case 2: this->z = val; break;
        case 3: this->w = val; break;
        default: return;
    }
}

void Point::wScale() {
    this->x = this->x / this->w;
    this->y = this->y / this->w;
    this->z = this->z / this->w;
    this->w = this->w / this->w;
}

void Point::print() const {
    stageFile << setprecision(7) << fixed << this->x << " " << this->y << " " << this->z << endl;
}

double Point::dotProduct(const Point& p) const {
    double a = this->x * p.x;
    double b = this->y * p.y;
    double c = this->z * p.z;
    return a + b + c;
}

Point Point::crossProduct(const Point& p) const {
    double a = y * p.z - z * p.y;
    double b = z * p.x - x * p.z;
    double c = x * p.y - y * p.x;
    return Point(a, b, c);
}

Point Point::scalarMultiply(double s) const {
    return Point(s * this->x, s * this->y, s * this->z);
}

void Point::normalize() {
    double length = sqrt( x*x + y*y + z*z );
    x = x / length;
    y = y / length;
    z = z / length;
}

Point Point::subtract(const Point& p) const {
    return Point(x - p.x, y - p.y, z - p.z);
}


class Triangle {
public:
    vector<Point> points;
    Color color;

    Triangle(vector<Point> points);
    double getMaxY();
    double getMinY();
    double getMaxX();
    double getMinX();
    vector<pair<double, double>> getColumns(double ys, double leftX, double dx);
};

Triangle::Triangle(vector<Point> points) {
    for (int i = 0; i < 3; i++) {
        this->points.push_back(points[i]);
    }
    color = Color(rand() % 256, rand() % 256, rand() % 256);
}

double Triangle::getMaxY() {
    return max(max(points[0].y, points[1].y), points[2].y);
}

double Triangle::getMinY() {
    return min(min(points[0].y, points[1].y), points[2].y);
}

double Triangle::getMaxX() {
    return max(max(points[0].x, points[1].x), points[2].x);
}

double Triangle::getMinX() {
    return min(min(points[0].x, points[1].x), points[2].x);
}

vector<pair<double, double>> Triangle::getColumns(double ys, double leftX, double dx) {
    double xs, zs;
    vector<pair<double, double>> columns, finalColumns;
    double minX = getMinX();
    double maxX = getMaxX();

    Point a, b, c;
    a = points[0];
    b = points[1];
    c = points[2];

    //    edge a - b
    double xa = a.x + ((ys - a.y)/(b.y - a.y))*(b.x - a.x);
    double za = a.z + ((ys - a.y)/(b.y - a.y))*(b.z - a.z);
    int s1 = (int) round(abs(leftX - xa) / dx);
    //    edge b - c
    double xb = b.x + ((ys - b.y)/(c.y - b.y))*(c.x - b.x);
    double zb = b.z + ((ys - b.y)/(c.y - b.y))*(c.z - b.z);
    int s2 = (int) round(abs(leftX - xb) / dx);
    //    edge c - a
    double xc = c.x + ((ys - c.y)/(a.y - c.y))*(a.x - c.x);
    double zc = c.z + ((ys - c.y)/(a.y - c.y))*(a.z - c.z);
    int s3 = (int) round(abs(leftX - xc) / dx);


    if (!isinf(xa) && !isinf(xb) && !isinf(xc)) {
        if (s1 == s2) {
            columns.push_back(make_pair(xa, za));
            columns.push_back(make_pair(xc, zc));
        }
        else if (s2 == s3) {
            columns.push_back(make_pair(xb, zb));
            columns.push_back(make_pair(xa, za));
        }
        else if (s3 == s1) {
            columns.push_back(make_pair(xc, zc));
            columns.push_back(make_pair(xb, zb));
        } else {
            columns.push_back(make_pair(xa, za));
            columns.push_back(make_pair(xb, zb));
            columns.push_back(make_pair(xc, zc));
        }
    } else {
        if (!isinf(xa)) columns.push_back(make_pair(xa, za));
        if (!isinf(xb)) columns.push_back(make_pair(xb, zb));
        if (!isinf(xc)) columns.push_back(make_pair(xc, zc));
    }

    for (int i = 0; i < columns.size(); i++) {
        xs = columns[i].first;
        zs = columns[i].second;
        if (xs >= minX && xs <= maxX) {
            finalColumns.push_back(make_pair(xs, zs));
        }
    }

    return finalColumns;
}


class TransformationMatrix {
public:
    double matrix[4][4];
    TransformationMatrix() {}
    TransformationMatrix(double matrix[4][4]);
    Point transformPoint(const Point& p);
    TransformationMatrix product(TransformationMatrix t);
    static TransformationMatrix identity();
    static TransformationMatrix translate(double tx, double ty, double tz);
    static TransformationMatrix scale(double sx, double sy, double sz);
    static TransformationMatrix rotate(double angle, double ax, double ay, double az);
};


TransformationMatrix::TransformationMatrix(double matrix[4][4]) {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            this->matrix[i][j] = matrix[i][j];
        }
    }
}

TransformationMatrix TransformationMatrix::identity() {
    TransformationMatrix t = TransformationMatrix();
    t.matrix[0][0] = 1, t.matrix[0][1] = 0, t.matrix[0][2] = 0, t.matrix[0][3] = 0;
    t.matrix[1][0] = 0, t.matrix[1][1] = 1, t.matrix[1][2] = 0, t.matrix[1][3] = 0;
    t.matrix[2][0] = 0, t.matrix[2][1] = 0, t.matrix[2][2] = 1, t.matrix[2][3] = 0;
    t.matrix[3][0] = 0, t.matrix[3][1] = 0, t.matrix[3][2] = 0, t.matrix[3][3] = 1;
    return t;
}

TransformationMatrix TransformationMatrix::translate(double tx, double ty, double tz) {
    TransformationMatrix t = TransformationMatrix();
    t.matrix[0][0] = 1, t.matrix[0][1] = 0, t.matrix[0][2] = 0, t.matrix[0][3] = tx;
    t.matrix[1][0] = 0, t.matrix[1][1] = 1, t.matrix[1][2] = 0, t.matrix[1][3] = ty;
    t.matrix[2][0] = 0, t.matrix[2][1] = 0, t.matrix[2][2] = 1, t.matrix[2][3] = tz;
    t.matrix[3][0] = 0, t.matrix[3][1] = 0, t.matrix[3][2] = 0, t.matrix[3][3] = 1;
    return t;
}

TransformationMatrix TransformationMatrix::scale(double sx, double sy, double sz) {
    TransformationMatrix t = TransformationMatrix();
    t.matrix[0][0] = sx, t.matrix[0][1] = 0, t.matrix[0][2] = 0, t.matrix[0][3] = 0;
    t.matrix[1][0] = 0, t.matrix[1][1] = sy, t.matrix[1][2] = 0, t.matrix[1][3] = 0;
    t.matrix[2][0] = 0, t.matrix[2][1] = 0, t.matrix[2][2] = sz, t.matrix[2][3] = 0;
    t.matrix[3][0] = 0, t.matrix[3][1] = 0, t.matrix[3][2] = 0, t.matrix[3][3] = 1;
    return t;
}

Point RodriguesFormula(const Point& rotateAxis, const Point& rotateVector, double angle) {
     Point p1 = rotateAxis.scalarMultiply(cos(angle));
     double temp = rotateVector.dotProduct(rotateAxis) * (1 - cos(angle));
     Point p2 = rotateVector.scalarMultiply(temp);
     Point p3 = rotateVector.crossProduct(rotateAxis).scalarMultiply(sin(angle));
     return Point(p1.x + p2.x + p3.x, p1.y + p2.y + p3.y, p1.z + p2.z + p3.z);
}

TransformationMatrix TransformationMatrix::rotate(double angle, double ax, double ay, double az) {
    Point a = Point(ax, ay, az);
    a.normalize();
    double angleRadian = angle * (PI / 180);

    Point c1 = RodriguesFormula(Point(1, 0, 0), a, angleRadian);
    Point c2 = RodriguesFormula(Point(0, 1, 0), a, angleRadian);
    Point c3 = RodriguesFormula(Point(0, 0, 1), a, angleRadian);

    TransformationMatrix t = TransformationMatrix();
    t.matrix[0][0] = c1.x, t.matrix[0][1] = c2.x, t.matrix[0][2] = c3.x, t.matrix[0][3] = 0;
    t.matrix[1][0] = c1.y, t.matrix[1][1] = c2.y, t.matrix[1][2] = c3.y, t.matrix[1][3] = 0;
    t.matrix[2][0] = c1.z, t.matrix[2][1] = c2.z, t.matrix[2][2] = c3.z, t.matrix[2][3] = 0;
    t.matrix[3][0] = 0, t.matrix[3][1] = 0, t.matrix[3][2] = 0, t.matrix[3][3] = 1;
    return t;
}

Point TransformationMatrix::transformPoint(const Point& p) {
    Point transformedPoint = Point();
    double a, b, sum;
    for (int i = 0; i < 4; i++) {
        sum = 0;
        for (int j = 0; j < 4; j++) {
            a = this->matrix[i][j];
            b = p.getCoordinate(j);
            sum += a * b;
        }
        transformedPoint.setCoordinate(i, sum);
    }
    transformedPoint.wScale();
    return transformedPoint;
}

TransformationMatrix TransformationMatrix::product(TransformationMatrix t) {
    TransformationMatrix res = TransformationMatrix();
    double a, b, sum;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            sum = 0;
            for (int k = 0; k < 4; k++) {
                a = this->matrix[i][k];
                b = t.matrix[k][j];
                sum += a * b;
            }
            res.matrix[i][j] = sum;
        }
    }

    return res;
}


vector<Point> loadScene() {
    ifstream sceneFile;
    sceneFile.open("./scene.txt");

    if (!sceneFile) {
        cerr << "Unable to open file scene.txt" << endl;
        exit(1);
    }

    stageFile.open("stage1.txt");
    if (!stageFile) {
        cerr << "Unable to open file stage1.txt" << endl;
        exit(1);
    }

    sceneFile >> eyeX >> eyeY >> eyeZ;
    sceneFile >> lookX >> lookY >> lookZ;
    sceneFile >> upX >> upY >> upZ;
    sceneFile >> fovY >> aspectRatio >> near >> far;

    vector<Point> points;
    stack<TransformationMatrix> matrixStack;
    matrixStack.push(TransformationMatrix::identity());
    stack<unsigned long> pushStack;
    string command;

    while (getline(sceneFile, command)) {
        stringstream ss(command);
        ss >> command;

        if (command == "triangle") {
            Point p1, p2, p3;
            sceneFile >> p1.x >> p1.y >> p1.z;
            sceneFile >> p2.x >> p2.y >> p2.z;
            sceneFile >> p3.x >> p3.y >> p3.z;

            TransformationMatrix t = matrixStack.top();
            Point t1 = t.transformPoint(p1);
            Point t2 = t.transformPoint(p2);
            Point t3 = t.transformPoint(p3);
            t1.print();
            t2.print();
            t3.print();
            stageFile << endl;

            points.push_back(t1);
            points.push_back(t2);
            points.push_back(t3);
        } else if (command == "translate") {
            double tx, ty, tz;
            sceneFile >> tx >> ty >> tz;
            TransformationMatrix t = TransformationMatrix::translate(tx, ty, tz);
            matrixStack.push(matrixStack.top().product(t));
        } else if (command == "scale") {
            double sx, sy, sz;
            sceneFile >> sx >> sy >> sz;
            TransformationMatrix s = TransformationMatrix::scale(sx, sy, sz);
            matrixStack.push(matrixStack.top().product(s));
        } else if (command == "rotate") {
            double angle, ax, ay, az;
            sceneFile >> angle >> ax >> ay >> az;
            TransformationMatrix r = TransformationMatrix::rotate(angle, ax, ay, az);
            matrixStack.push(matrixStack.top().product(r));
        } else if (command == "push") {
            pushStack.push(matrixStack.size());
        } else if (command == "pop") {
            int i = 0;
            unsigned long popCount = matrixStack.size() - pushStack.top();
            pushStack.pop();

            while (i < popCount) {
                matrixStack.pop();
                i++;
            }
        } else if (command == "end") {
            break;
        } else {
            continue;
        }
    }

    sceneFile.close();
    stageFile.close();
    return points;
}

vector<Point> viewTransformation(const vector<Point>& points) {
    Point look = Point(lookX, lookY, lookZ);
    Point eye = Point(eyeX, eyeY, eyeZ);
    Point up = Point(upX, upY, upZ);
    Point l = look.subtract(eye);
    l.normalize();

    Point r = l.crossProduct(up);
    r.normalize();

    Point u = r.crossProduct(l);
    u.normalize();

    TransformationMatrix translation = TransformationMatrix().translate(-eyeX, -eyeY, -eyeZ);
    double temp[4][4] = {
            {r.x, r.y, r.z, 0},
            {u.x, u.y, u.z, 0},
            {-l.x, -l.y, -l.z, 0},
            {0, 0, 0, 1}
    };
    TransformationMatrix rotation = TransformationMatrix(temp);
    TransformationMatrix v = rotation.product(translation);

    stageFile.open("stage2.txt");
    if (!stageFile) {
        cerr << "Unable to open file stage2.txt" << endl;
        exit(1);
    }
    vector<Point> transformedPoints;
    for (int i = 0; i < points.size(); i += 3) {
        for (int j = 0; j < 3; j++) {
            Point p = v.transformPoint(points[i+j]);
            p.print();
            transformedPoints.push_back(p);
        }
        stageFile << endl;
    }
    stageFile.close();
    return transformedPoints;
}

vector<Point> projectionTransformation(vector<Point> points) {
    double fovX = fovY * aspectRatio;
    double t = near * tan((fovY / 2)  * (PI / 180));
    double r = near * tan((fovX / 2)  * (PI / 180));

    double temp[4][4] = {
            {near/r, 0, 0, 0},
            {0, near/t, 0, 0},
            {0, 0, -(far+near)/(far-near), -(2*far*near)/(far-near)},
            {0, 0, -1, 0}
    };
    TransformationMatrix projection = TransformationMatrix(temp);

    stageFile.open("stage3.txt");
    if (!stageFile) {
        cerr << "Unable to open file stage3.txt" << endl;
        exit(1);
    }
    vector<Point> transformedPoints;
    for (int i = 0; i < points.size(); i += 3) {
        for (int j = 0; j < 3; j++) {
            Point p = projection.transformPoint(points[i+j]);
            p.print();
            transformedPoints.push_back(p);
        }
        stageFile << endl;
    }
    stageFile.close();
    return transformedPoints;
}

void loadConfig() {
    ifstream configFile;
    configFile.open("./config.txt");

    if (!configFile) {
        cerr << "Unable to open file config.txt" << endl;
        exit(1);
    }

    configFile >> screenWidth >> screenHeight;
    configFile >> xLimit;
    configFile >> yLimit;
    configFile >> zFrontLimit >> zRearLimit;
    configFile.close();
}

vector<Triangle> generateTriangles(vector<Point> points) {

    vector<Triangle> triangles;

    for (int i = 0; i < points.size(); i += 3) {
        vector<Point> triPoints;
        for (int j = 0; j < 3; j++) {
            triPoints.push_back(points[i+j]);
        }
        Triangle t = Triangle(triPoints);
        triangles.push_back(t);
    }
    return triangles;
}

void zBufferAlgorithm(vector<Point> points) {

    vector<Triangle> triangles = generateTriangles(points);


    double xLeft = xLimit;
    double xRight = -xLimit;
    double dx = (xRight - xLeft) / screenWidth;

    double yBottom = yLimit;
    double yTop = -yLimit;
    double dy = (yTop - yBottom) / screenHeight;

    double topY = yTop - dy / 2;
    double bottomY = yBottom + dy / 2;
    double leftX = xLeft + dx / 2;
    double rightX = xRight - dx / 2;

    vector<vector<double>> zBuffer;
    vector<double> zRowVector;
    bitmap_image image(screenWidth, screenHeight);

    for (int i = 0; i < screenHeight; i++) {
        zRowVector.clear();
        for (int j = 0; j < screenWidth; j++) {
            zRowVector.push_back(zRearLimit);
            image.set_pixel(i, j, 0, 0, 0);
        }
        zBuffer.push_back(zRowVector);
    }

    int top_scanline, bottom_scanline, left_scanline, right_scanline;
    double maxY, minY, maxX, minX;
    for (int i = 0; i < triangles.size(); i++) {
        Triangle triangle = triangles[i];

        maxY = triangle.getMaxY();
        if (maxY > topY) {
            top_scanline = 0;
        } else {
            top_scanline = (int) round(abs(topY - maxY) / dy);
        }

        minY = triangle.getMinY();
        if (minY < bottomY) {
            bottom_scanline = screenHeight - 1;
        } else {
            bottom_scanline = (int) round(abs(topY - minY) / dy);
        }

        cout << "Triangle " << i << endl;
        cout << "top - " << top_scanline << "\tbottom - " << bottom_scanline << endl;

        for (int row = top_scanline; row <= bottom_scanline; row++) {
            double ys = topY - row * dy;
            vector<pair<double, double>> columns = triangle.getColumns(ys, leftX, dx);

            double za, zb;
            if (columns.empty()) continue;
            if (columns.size() == 1) {
                minX = columns[0].first;
                maxX = columns[0].first;
                za = columns[0].second;
                zb = columns[0].second;
            } else {
                if (columns[0].first <= columns[1].first) {
                    minX = columns[0].first;
                    maxX = columns[1].first;
                    za = columns[0].second;
                    zb = columns[1].second;
                } else {
                    minX = columns[1].first;
                    maxX = columns[0].first;
                    za = columns[1].second;
                    zb = columns[0].second;
                }
            }

            if (minX < leftX) {
                left_scanline = 0;
            } else {
                left_scanline = (int) round(abs(leftX - minX) / dx);
            }

            if (maxX > rightX) {
                right_scanline = screenWidth - 1;
            } else {
                right_scanline = (int) round(abs(leftX - maxX) / dx);
            }
            cout << "left - " << left_scanline << "\tright - " << right_scanline << endl;

            double z, inc;
            z = za;
            if (right_scanline > left_scanline) {
                inc = (zb - za) / (right_scanline - left_scanline);
            } else {
                inc = 0;
            }

            if (maxX - minX != 0) {
                for (int col = left_scanline; col <= right_scanline; col++) {
                    if (z >= zFrontLimit && z <= zRearLimit && z < zBuffer[row][col]) {
                        zBuffer[row][col] = z;
                        Color c = triangle.color;
                        image.set_pixel(col, row, c.r, c.g, c.b);
                    }
                    z += inc;
                }
            }
        }
    }

    // save the image and z buffer
    image.save_image("output.bmp");
    ofstream bufferFile;
    bufferFile.open("z_buffer.txt");
    for (int row = 0; row < zBuffer.size(); row++) {
        for (int col = 0; col < zBuffer[row].size(); col++) {
            if (zBuffer[row][col] < zRearLimit) {
                bufferFile << zBuffer[row][col] << "\t";
            }
        }
        bufferFile << endl;
    }

    // clear the buffer and the image
    for (int i = 0; i < zBuffer.size(); i++) {
        zBuffer[i].clear();
    }
    zBuffer.clear();
    image.clear();
}


int main() {

    vector<Point> stage1Points = loadScene();
    vector<Point> stage2Points = viewTransformation(stage1Points);
    vector<Point> stage3Points = projectionTransformation(stage2Points);

    loadConfig();
    zBufferAlgorithm(stage3Points);
    return 0;
}