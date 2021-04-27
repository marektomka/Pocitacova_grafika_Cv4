#pragma once
#include <QtWidgets>

class W_edge;
class H_edge;

class Vertex {
private:
	double x, y, z;
	H_edge* edge;
public:
	Vertex() { x = 0, y = 0, z = 0; }
	Vertex(double x1, double y1, double z1) { x = x1, y = y1, z = z1; }
	void setCoor(double x1, double y1, double z1) { x = x1; y = y1; z = z1; }
	double X() { return x; }
	double Y() { return y; }
	double Z() { return z; }

};

class Face {
	H_edge* edge;
public:
	Face() { edge = nullptr; }
	void setPointer(H_edge* edge1) { edge = edge1; }

};
/*
class Face {
	W_edge* edge;
	int side;
public:
	Face() { edge = nullptr; }
	void setEdge(W_edge* edge1, int sd) { edge = edge1, side = sd; }

};*/


class H_edge {
	Vertex* vert_origin;
	Face* face;
	H_edge* edge_prev, * edge_next;
	H_edge* pair;
public:
	H_edge() { vert_origin = nullptr, face = nullptr, edge_prev = nullptr, edge_next = nullptr, pair = nullptr; }
	void setAllPtr(Vertex* origin, Face* fac, H_edge* prev, H_edge* next, H_edge* parik) { vert_origin = origin; face = fac; edge_prev = prev; edge_next = next; pair = parik; }
};

/*
class W_edge {
	Vertex* vert_origin, * vert_destination;
	Face* face_left, * face_right;
	W_edge* edge_left_prev, * edge_left_next;
	W_edge* edge_right_prev, * edge_right_next;
public:
	W_edge() { vert_origin = nullptr, vert_destination = nullptr, face_left = nullptr, face_right = nullptr, edge_left_prev = nullptr, edge_left_next = nullptr, edge_right_prev = nullptr, edge_right_next = nullptr; }
	void setVertexPtr(Vertex* origin, Vertex* destin) { vert_origin = origin; vert_destination = destin; }
	void setFacePtr(Face* fleft, Face* fright) { face_left = fleft; face_right = fright; }
	void setLeftEdgePtr(W_edge* leftPrev, W_edge* leftNext) { edge_left_prev = leftPrev; edge_left_next = leftNext; }
	void setRightEdgePtr(W_edge* rightPrev, W_edge* rightNext) { edge_right_prev = rightPrev; edge_right_next = rightNext; }
};*/


class ViewerWidget :public QWidget {
	Q_OBJECT
private:
	QString name = "";
	QSize areaSize = QSize(0, 0);
	QImage* img = nullptr;
	QRgb* data = nullptr;
	QPainter* painter = nullptr;
	QVector<QPoint> points;
	QVector<QPoint> end; 
	QVector<QPoint> newend;
	QVector<QPoint> Pd;
	double dt = 0,degree = 0;
	int cpoint = 0,method = 0;
	QColor color;
	double F0 = 0, F1 = 0, F2 = 0, F3 = 0;
	double B0 = 0, B1 = 0, B2 = 0, B3 = 0;

	//---------------------------------------- Sphere var--------------------------------------//
	QVector<Vertex> tableVertex;
	QVector<H_edge> Hedges = QVector<H_edge>();
	QVector<Face> faces = QVector<Face>();
	//QVector<QVector<Face>> tableFace;
	//QVector<QVector<W_edge>> tableWEdge;
	double Fi = (1.0 + sqrt(5)) / 2;

public:
	ViewerWidget(QString viewerName, QSize imgSize, QWidget* parent = Q_NULLPTR);
	~ViewerWidget();
	void resizeWidget(QSize size);

	//Image functions
	bool setImage(const QImage& inputImg);
	QImage* getImage() { return img; };
	bool isEmpty();

	//Data functions
	QRgb* getData() { return data; }
	void setPixel(int x, int y, const QColor& color);
	void setPixel(int x, int y, unsigned char r, unsigned char g, unsigned char b);
	bool isInside(int x, int y) { return (x >= 0 && y >= 0 && x < img->width() && y < img->height()) ? true : false; }

	//Get/Set functions
	QString getName() { return name; }
	void setName(QString newName) { name = newName; }

	void setPainter() { painter = new QPainter(img); }
	void setDataPtr() { data = reinterpret_cast<QRgb*>(img->bits()); }

	int getImgWidth() { return img->width(); };
	int getImgHeight() { return img->height(); };

	void clear(QColor color = Qt::white);

	void drawDDA(int x1, int y1, int x2, int y2,QColor col);
	void drawBresen(int x1, int y1, int x2, int y2, QColor col);
	void drawBresenCircle(int x1, int y1, int x2, int y2, QColor col);
	void drawCurve(QColor col, int meth);
	void setPoints(QPoint p);
	void clearButton();
	void drawHermite();
	void drawCoons();
	void drawBezier();
	void setFpolyn(double te);
	void setBpolyn(double te);
	void setHermiteFactor(double d,int p);
	void drawPoints();
	
	//------------------------------------Sphere functions-------------------------------------//
	void createSphere();
	void fillVertex();
	void fillHedges();
	void fillFaces();

public slots:
	void paintEvent(QPaintEvent* event) Q_DECL_OVERRIDE;
};