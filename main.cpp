#include <bits/stdc++.h>

using namespace std;

#define PI 3.14159265358979323846

double eyeX, eyeY, eyeZ;
double lookX, lookY, lookZ;
double upX, upY, upZ;
double fovY, aspectRatio, near, far;

double screenWidth, screenHeight;
double xLimit, yLimit;
double zFrontLimit, zRearLimit;

ofstream stageFile;

class Point {
public:
    double x{}, y{}, z{}, w;
    Point();
    Point(double x, double y, double z);
    double getCoordinate(int index) const;
    void setCoordinate(int index, double val);
    void wScale();
    void print() const;
    double dotProduct(Point p) const;
    Point crossProduct(Point p) const;
    Point scalarMultiply(double s) const;
    Point subtract(Point p) const;
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

double Point::dotProduct(Point p) const {
    double a = this->x * p.x;
    double b = this->y * p.y;
    double c = this->z * p.z;
    return a + b + c;
}

Point Point::crossProduct(Point p) const {
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

Point Point::subtract(Point p) const {
    return Point(x - p.x, y - p.y, z - p.z);
}


class Triangle {
public:
    vector<Point> points;
    vector<int> color;

    Triangle(vector<Point> points);
};

Triangle::Triangle(vector<Point> points) {
    for (int i = 0; i < 3; i++) {
        this->points.push_back(points[i]);
    }
    for (int i = 0; i < 3; i++) {
        this->color.push_back(rand() % 256);
    }
}


class TransformationMatrix {
public:
    double matrix[4][4];
    TransformationMatrix() {}
    TransformationMatrix(double matrix[4][4]);
    Point transformPoint(Point p);
    TransformationMatrix product(TransformationMatrix t);
    void print();
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

Point RodriguesFormula(Point rotateAxis, Point rotateVector, double angle) {
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

Point TransformationMatrix::transformPoint(Point p) {
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

void TransformationMatrix::print() {
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << this->matrix[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
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
}

vector<Triangle> generateTriangles(vector<Point> points) {
    vector<Triangle> triangles;

    for (int i = 0; i < points.size(); i += 3) {
        vector<Point> triPoints(3);
        for (int j = 0; j < 3; j++) {
            triPoints.push_back(points[i+j]);
        }
        Triangle t = Triangle(triPoints);
        triangles.push_back(t);
    }
    return triangles;
}


int main() {

    vector<Point> stage1Points = loadScene();
    vector<Point> stage2Points = viewTransformation(stage1Points);
    vector<Point> stage3Points = projectionTransformation(stage2Points);

    loadConfig();
    vector<Triangle> triangles = generateTriangles(stage3Points);
    return 0;
}