#include   "ViewerWidget.h"
#include <math.h> 

ViewerWidget::ViewerWidget(QString viewerName, QSize imgSize, QWidget* parent)
	: QWidget(parent)
{
	setAttribute(Qt::WA_StaticContents);
	setMouseTracking(true);
	name = viewerName;
	if (imgSize != QSize(0, 0)) {
		img = new QImage(imgSize, QImage::Format_ARGB32);
		img->fill(Qt::white);
		resizeWidget(img->size());
		setPainter();
		setDataPtr();
	}
}
ViewerWidget::~ViewerWidget()
{
	delete painter;
	delete img;
}
void ViewerWidget::resizeWidget(QSize size)
{
	this->resize(size);
	this->setMinimumSize(size);
	this->setMaximumSize(size);
}

//Image functions
bool ViewerWidget::setImage(const QImage& inputImg)
{
	if (img != nullptr) {
		delete img;
	}
	img = new QImage(inputImg);
	if (!img) {
		return false;
	}
	resizeWidget(img->size());
	setPainter();
	update();

	return true;
}
bool ViewerWidget::isEmpty()
{
	if (img->size() == QSize(0, 0)) {
		return true;
	}
	return false;
}

//Data function
void ViewerWidget::setPixel(int x, int y, const QColor& color)
{
	if (isInside(x, y)) {
		data[x + y * img->width()] = color.rgb();
	}
}
void ViewerWidget::setPixel(int x, int y, unsigned char r, unsigned char g, unsigned char b)
{
	if (isInside(x, y)) {
		QColor color(r, g, b);
		setPixel(x, y, color);
	}
}

void ViewerWidget::clear(QColor color)
{
	for (size_t x = 0; x < img->width(); x++)
	{
		for (size_t y = 0; y < img->height(); y++)
		{
			setPixel(x, y, color);
		}
	}
	update();
}


void ViewerWidget::drawDDA(int x1, int y1, int x2, int y2, QColor col)
{	


	double m = DBL_MAX;      // smernica
	int X, Y;               // pomocne pre vymenu suradnic
	double xo, yo;			// pomocne pre posunutie v danom smere


	m = (y2 - y1) / (double)(x2 - x1);
	m = abs(m);

	if (m <= 1)
	{
		if (x1 > x2)
		{
			X = x2; x2 = x1; x1 = X;
			Y = y2; y2 = y1; y1 = Y;
		}
		yo = y1;
		setPixel(x1, y1, col);
		
		do
		{
			x1 = x1 + 1;
			if (y1 < y2) {
				yo = yo + m;
			}
			else {
				yo = yo - m;
			}
			setPixel(x1, static_cast<int>(yo + 0.5), col);
		} while (x1 < x2);
	}
	else if (m >= 1) 
	{
		if (y1 > y2)
		{
			Y = y2; y2 = y1; y1 = Y;
			X = x2; x2 = x1; x1 = X;
		}
		xo = x1;
		setPixel(x1, y1, col);

		do
		{
			y1 = y1 + 1;
			if (x1 < x2) {
				xo = xo + (1 / m);
			}
			else {
				xo = xo - (1 / m);
			}
			setPixel(static_cast<int>(xo + 0.5), y1, col);
		} while (y1 < y2);
	}
	update();
}

void ViewerWidget::drawBresen(int x1, int y1, int x2, int y2, QColor col)
{
	int k1, k2, p1;
	int xo, yo;
	double m = DBL_MAX;
	int X, Y;


	m = (y2 - y1) / (double)(x2 - x1);

	if (m > 0 && m < 1) 
	{
		if (x1 > x2) 
		{
			X = x2; x2 = x1; x1 = X;
			Y = y2; y2 = y1; y1 = Y;
		}
		k1 = 2 * (y2 - y1);
		k2 = k1 - 2 * (x2 - x1);
		p1 = k1 - (x2 - x1);
		setPixel(x1, y1, col);

		do
		{
			x1 += 1;
			if (p1 > 0) {
			p1 = p1 + k2;
				if (y1 < y2) {
					y1 += 1;
				}
				else {
					y1 -= 1;
				}
			}
			else {
				p1 = p1 + k1;
			}
			setPixel(x1, y1, col);
		} while (x1 < x2);
	}
	else if (m > -1 && m <= 0) 
	{
		if (x1 > x2)
		{
			X = x2; x2 = x1; x1 = X;
			Y = y2; y2 = y1; y1 = Y;
		}
		k1 = 2 * (y2 - y1);
		k2 = k1 + 2 * (x2 - x1);
		p1 = k1 + (x2 - x1);
		setPixel(x1, y1, col);

		do
		{
			x1 += 1;
			if (p1 < 0) {
				y1 -= 1;
				p1 = p1 + k2;
			}
			else {
				p1 = p1 + k1;
			}
			setPixel(x1, y1, col);

		} while (x1 < x2);

	}
	else if (m >= 1) 
	{
		if (y1 > y2)
		{
			X = x2; x2 = x1; x1 = X;
			Y = y2; y2 = y1; y1 = Y;
		}
		k1 = 2 * (x2 - x1);
		k2 = k1 - 2 * (y2 - y1);
		p1 = k1 - (y2 - y1);
		setPixel(x1, y1, col);

		do
		{
			y1 += 1;
			if (p1 > 0) {
				p1 = p1 + k2;
				if (y1 < y2) {
					x1 += 1;
				}
				else {
					x1 -= 1;
				}
			}
			else {
				p1 = p1 + k1;
			}
			setPixel(x1, y1, col);

		} while (y1 < y2);

	}
	else if (m <= -1)
	{
		if (y1 > y2)
		{
			X = x2; x2 = x1; x1 = X;
			Y = y2; y2 = y1; y1 = Y;
		}

		k1 = 2 * (x2 - x1);
		k2 = k1 + 2 * (y2 - y1);
		p1 = k1 + (y2 - y1);
		setPixel(x1, y1, col);

		do
		{
			y1 += 1;
			if (p1 < 0) {
				x1 -= 1;
				p1 = p1 + k2;
			}
			else {
				p1 = p1 + k1;
			}
			setPixel(x1, y1, col);
		} while (y1 < y2);
	}
	update();
}

void ViewerWidget::drawBresenCircle(int x1, int y1, int x2, int y2, QColor col)
{
	int r = 0;
	int p1;
	int X, Y,dvaX,dvaY;
	r = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
	p1 = 1 - r;
	X = 0;
	Y = r ;
	dvaX = 3;
	dvaY = 2 * r + 2;
	do
	{
		setPixel(X+x1, Y+y1, col);
		setPixel(-X + x1, Y + y1, col);
		setPixel(X + x1, -Y + y1, col);
		setPixel(-X + x1, -Y + y1, col);
		setPixel(Y + x1,X + y1, col);
		setPixel(-Y + x1, X + y1, col);
		setPixel(Y + x1, -X + y1, col);
		setPixel(-Y + x1, -X + y1, col);
		if (p1 > 0) {
			p1 = p1 - dvaY;
			Y -= 1;
			dvaY -= 2;
		}
		p1 = p1 + dvaX;
		dvaX += 2;
		X += 1;
		

	} while (X<=Y);

	update();
}

void ViewerWidget::drawCurve(QColor col, int meth)
{
	clear();
	color = col;
	method = meth;
	if (method == 0) {
		if (points.size() >= 2) {
			drawHermite();
		}
	}
	else if (method == 1) {
		if (points.size() >= 3) {
			drawBezier();
		}
	}
	else if (method == 2) {
		if (points.size() >= 4) {
			drawCoons();
		}
	}

}

void ViewerWidget::setPoints(QPoint p)
{
	points.push_back(p);
}

void ViewerWidget::clearButton()
{
	points.clear();
	newend.clear();
	end.clear();
	degree = 0;
	cpoint = 0;
	dt = 0;
	F0 = 0, F1 = 0, F2 = 0, F3 = 0;
}

void ViewerWidget::drawHermite()
{
	QVector<QPoint> P;		// P - body .... Pd - dotycnice.... end - koncove body dotycnic
	QPointF Q0, Q1;
	double t = 0;
	P = points;
	dt = 0.05;

	if (degree == 0 && cpoint == 0) {
		for (int i = 0; i < P.size(); i++)									// vypocet dotycnic
		{
			end.append(QPoint(P[i].x(), P[i].y() + 100));
			Pd.append(QPoint(end[i] - P[i]));
			newend.push_back(end[i]);
		}

	}
	// s tymto neviem nejak zabezpecit, aby som mohol menit viacej dotycnic pri jednej krivke... 
	//end[cpoint].setX(P[cpoint].x() + 100*sin(degree * (M_PI / 180)));		// vypocet dotycnic so stupnami 
	//end[cpoint].setY(P[cpoint].y() + 100*cos(degree * (M_PI / 180)));

	end[cpoint].setX(((newend[cpoint].x() - P[cpoint].x()) * cos(degree * (M_PI / 180))) - ((newend[cpoint].y() - P[cpoint].y()) * sin(degree * (M_PI / 180))) + P[cpoint].x());      // pouzite z rotacie usecky z predosleho zadania
	end[cpoint].setY(((newend[cpoint].x() - P[cpoint].x()) * sin(degree * (M_PI / 180))) + ((newend[cpoint].y() - P[cpoint].y()) * cos(degree * (M_PI / 180))) + P[cpoint].y());

	Pd[cpoint] = end[cpoint] - P[cpoint];

	for (int i = 0; i < P.size(); i++)
	{
		drawDDA(P[i].x(), P[i].y(), end[i].x(),end[i].y(), Qt::red);
	}
	
	for (double i = 1; i < P.size(); i++)
	{
		Q0 = P[i - 1];
		t = dt;
		do{
			setFpolyn(t);
			Q1 = P[i - 1] * F0 + P[i] * F1 + Pd[i - 1] * F2 + Pd[i] * F3;
			drawDDA(int(Q0.x() + 0.5), int(Q0.y() + 0.5), int(Q1.x() + 0.5), int(Q1.y() + 0.5),color);
			Q0 = Q1;
			t += dt;
		} while (t < 1);

		drawDDA(int(Q0.x() + 0.5), int(Q0.y() + 0.5), int(P[i].x() + 0.5), int(P[i].y() + 0.5), color);
	}
	
	drawPoints();

}

void ViewerWidget::drawCoons()
{
	QVector<QPoint> P = points;
	QPointF Q0, Q1;
	double t;
	dt = 0.05;
	for  (int i = 3; i < P.size(); i++)
	{
		t = 0;
		setBpolyn(t);
		Q0 = P[i - 3] * B0 + P[i - 2] * B1 + P[i - 1] * B2 + P[i] * B3;
		
		do{
			t += dt;
			setBpolyn(t);
			Q1 = P[i - 3] * B0 + P[i - 2] * B1 + P[i - 1] * B2 + P[i] * B3;
			drawDDA(int(Q0.x()+0.5), int(Q0.y()+0.5), int(Q1.x()+0.5), int(Q1.y()+0.5), color);
			Q0 = Q1;
		} while (t < 1);

	}
	drawPoints();
}

void ViewerWidget::drawBezier()
{
	QVector<QPoint> P = points;
	QPointF Q0, Q1;
	QVector<QVector<QPoint>> field;
	field.resize(P.size());
	for (int i = 0; i < P.size(); i++){
		field[0].append(P[i]);
	}
	double t;
	dt = 0.05;
	t = dt;
	Q0 = P[0];

	do
	{
		for (int i = 1; i < P.size(); i++){
			for (int j = 0; j < P.size() - i ; j ++){
				field[i].append(QPoint((1.0 - t) * field[i - 1][j] + t * field[i - 1][j + 1]));
			}
		}
		Q1 = field[P.size() - 1][0];
		drawDDA(int(Q0.x()+0.5), int(Q0.y() + 0.5), int(Q1.x() + 0.5), int(Q1.y() + 0.5), color);
		Q0 = Q1;
		t += dt;

		field.clear();
		field.resize(P.size());
		for (int i = 0; i < P.size(); i++){
			field[0].append(P[i]);
		}
		
	} while (t < 1);

	drawDDA(int(Q0.x() + 0.5), int(Q0.y() + 0.5), P[P.size() - 1].x(), P[P.size() - 1].y(), color);
	drawPoints();
}

void ViewerWidget::setFpolyn(double te)
{
	F0 = 2 * pow(te, 3) - 3 * pow(te, 2) + 1;
	F1 = -2 * pow(te, 3) + 3 * pow(te, 2);
	F2 = pow(te, 3) - 2 * pow(te, 2) + te;
	F3 = pow(te, 3) - pow(te, 2);
}

void ViewerWidget::setBpolyn(double te)
{
	
	B0 = (-1 / 6.0) * pow(te, 3) + (1 / 2.0) * pow(te, 2) - (1 / 2.0) * te + (1 / 6.0);
	B1 = (1 / 2.0) * pow(te, 3) - pow(te, 2) + (2 / 3.0);
	B2 = (-1 / 2.0) * pow(te, 3) + (1 / 2.0) * pow(te, 2) + (1 / 2.0) * te + (1 / 6.0);
	B3 = (1 / 6.0) * pow(te, 3);
}

void ViewerWidget::setHermiteFactor(double d,int p)
{
	degree = d;
	cpoint = p;
	drawCurve(color,method);
}

void ViewerWidget::drawPoints()
{
	for (int i = 0; i < points.size(); i++)
	{
		setPixel(points[i].x(), points[i].y(), Qt::blue);
		setPixel(points[i].x() + 1, points[i].y(), Qt::blue);
		setPixel(points[i].x(), points[i].y() + 1, Qt::blue);
		setPixel(points[i].x() - 1, points[i].y(), Qt::blue);
		setPixel(points[i].x(), points[i].y() - 1 , Qt::blue);
	}
	update();
}

//---------------------------------------- Sphere stuff-----------------------------------------------//

void ViewerWidget::createSphere()
{
	Hedges.resize(60);
	faces.resize(20);

	fillVertex();
	fillFaces();
	fillHedges();
	
	/*for (int i = 0; i < faces.size(); i++)
	{
		faces[i].setPointer(&Hedges[i]);
	}
	for (size_t i = 0; i < tableVertex.size(); i++)
	{
		qDebug() << tableVertex[i].X() << " " << tableVertex[i].Y() << " " << tableVertex[i].Z();
	}
	*/
	
}

void ViewerWidget::fillVertex()
{
	// iba suradnice zatial // 
	tableVertex.append(Vertex(0.0, 1.0 / sqrt(1. + Fi * Fi), Fi / sqrt(1.0 + Fi * Fi)));
	tableVertex.append(Vertex(0.0, 1.0 / sqrt(1. + Fi * Fi), -Fi / sqrt(1.0 + Fi * Fi)));
	tableVertex.append(Vertex(0.0, -1.0 / sqrt(1. + Fi * Fi), Fi / sqrt(1.0 + Fi * Fi)));
	tableVertex.append(Vertex(0.0, -1.0 / sqrt(1. + Fi * Fi), -Fi / sqrt(1.0 + Fi * Fi)));

	tableVertex.append(Vertex(1.0 / sqrt(1.0 + Fi * Fi), Fi / sqrt(1.0 + Fi * Fi), 0.0));
	tableVertex.append(Vertex(1.0 / sqrt(1.0 + Fi * Fi), -Fi / sqrt(1.0 + Fi * Fi), 0.0));
	tableVertex.append(Vertex(-1.0 / sqrt(1.0 + Fi * Fi), Fi / sqrt(1.0 + Fi * Fi), 0.0));
	tableVertex.append(Vertex(-1.0 / sqrt(1.0 + Fi * Fi), -Fi / sqrt(1.0 + Fi * Fi), 0.0));

	tableVertex.append(Vertex(Fi / sqrt(1.0 + Fi * Fi), 0.0, 1.0 / sqrt(1.0 + Fi * Fi)));
	tableVertex.append(Vertex(Fi / sqrt(1.0 + Fi * Fi), 0.0, -1.0 / sqrt(1.0 + Fi * Fi)));
	tableVertex.append(Vertex(-Fi / sqrt(1.0 + Fi * Fi), 0.0, 1.0 / sqrt(1.0 + Fi * Fi)));
	tableVertex.append(Vertex(-Fi / sqrt(1.0 + Fi * Fi), 0.0, -1.0 / sqrt(1.0 + Fi * Fi)));

}

void ViewerWidget::fillHedges()
{
	
	Hedges[0].setAllPtr(&tableVertex[1],&faces[0],&Hedges[4], &Hedges[2], &Hedges[1]);
	Hedges[1].setAllPtr(&tableVertex[0], &faces[1], &Hedges[8], &Hedges[6], &Hedges[0]);
	Hedges[2].setAllPtr(&tableVertex[0], &faces[0], &Hedges[0], &Hedges[4], &Hedges[3]);
	Hedges[3].setAllPtr(&tableVertex[4], &faces[7], &Hedges[36], &Hedges[35], &Hedges[2]);
	Hedges[4].setAllPtr(&tableVertex[4], &faces[0], &Hedges[2], &Hedges[0], &Hedges[5]);
	Hedges[5].setAllPtr(&tableVertex[1], &faces[4], &Hedges[15], &Hedges[45], &Hedges[4]);
	Hedges[6].setAllPtr(&tableVertex[0], &faces[1], &Hedges[1], &Hedges[8], &Hedges[7]);
	Hedges[7].setAllPtr(&tableVertex[6], &faces[2], &Hedges[12], &Hedges[10], &Hedges[6]);
	Hedges[8].setAllPtr(&tableVertex[6], &faces[1], &Hedges[6], &Hedges[1], &Hedges[9]);
	Hedges[9].setAllPtr(&tableVertex[0], &faces[9], &Hedges[30], &Hedges[29], &Hedges[8]);
	Hedges[10].setAllPtr(&tableVertex[1], &faces[2], &Hedges[7], &Hedges[12], &Hedges[11]);
	Hedges[11].setAllPtr(&tableVertex[11], &faces[3], &Hedges[16], &Hedges[14], &Hedges[10]);
	Hedges[12].setAllPtr(&tableVertex[11], &faces[2], &Hedges[10], &Hedges[7], &Hedges[13]);
	Hedges[13].setAllPtr(&tableVertex[6], &faces[11], &Hedges[24], &Hedges[23], &Hedges[12]);
	Hedges[14].setAllPtr(&tableVertex[1], &faces[3], &Hedges[11], &Hedges[16], &Hedges[15]);
	Hedges[15].setAllPtr(&tableVertex[9], &faces[4], &Hedges[45], &Hedges[5], &Hedges[14]);
	Hedges[16].setAllPtr(&tableVertex[9], &faces[3], &Hedges[14], &Hedges[11], &Hedges[17]);
	Hedges[17].setAllPtr(&tableVertex[11], &faces[13], &Hedges[20], &Hedges[18], &Hedges[16]);
	Hedges[18].setAllPtr(&tableVertex[9], &faces[13], &Hedges[17], &Hedges[20], &Hedges[19]);
	Hedges[19].setAllPtr(&tableVertex[3], &faces[14], &Hedges[46], &Hedges[43], &Hedges[18]);
	Hedges[20].setAllPtr(&tableVertex[3], &faces[13], &Hedges[18], &Hedges[17], &Hedges[21]);
	Hedges[21].setAllPtr(&tableVertex[11], &faces[12], &Hedges[22], &Hedges[59], &Hedges[20]);
	Hedges[22].setAllPtr(&tableVertex[7], &faces[12], &Hedges[59], &Hedges[21], &Hedges[23]);
	Hedges[23].setAllPtr(&tableVertex[11], &faces[11], &Hedges[13], &Hedges[24], &Hedges[22]);
	Hedges[24].setAllPtr(&tableVertex[7], &faces[11], &Hedges[23], &Hedges[13], &Hedges[25]);
	Hedges[25].setAllPtr(&tableVertex[6], &faces[10], &Hedges[28], &Hedges[26], &Hedges[24]);
	Hedges[26].setAllPtr(&tableVertex[7], &faces[10], &Hedges[25], &Hedges[28], &Hedges[27]);
	Hedges[27].setAllPtr(&tableVertex[10], &faces[18], &Hedges[55], &Hedges[56], &Hedges[26]);
	Hedges[28].setAllPtr(&tableVertex[10], &faces[10], &Hedges[26], &Hedges[25], &Hedges[29]);
	Hedges[29].setAllPtr(&tableVertex[6], &faces[9], &Hedges[9], &Hedges[30], &Hedges[28]);
	Hedges[30].setAllPtr(&tableVertex[10], &faces[9], &Hedges[29], &Hedges[9], &Hedges[31]);
	Hedges[31].setAllPtr(&tableVertex[0], &faces[8], &Hedges[34], &Hedges[32], &Hedges[30]);
	Hedges[32].setAllPtr(&tableVertex[10], &faces[8], &Hedges[31], &Hedges[34], &Hedges[33]);
	Hedges[33].setAllPtr(&tableVertex[8], &faces[17], &Hedges[53], &Hedges[54], &Hedges[32]);
	Hedges[34].setAllPtr(&tableVertex[8], &faces[8], &Hedges[32], &Hedges[31], &Hedges[35]);
	Hedges[35].setAllPtr(&tableVertex[0], &faces[7], &Hedges[3], &Hedges[36], &Hedges[34]);
	Hedges[36].setAllPtr(&tableVertex[8], &faces[7], &Hedges[35], &Hedges[3], &Hedges[37]);
	Hedges[37].setAllPtr(&tableVertex[4], &faces[6], &Hedges[40], &Hedges[38], &Hedges[36]);
	Hedges[38].setAllPtr(&tableVertex[8], &faces[6], &Hedges[37], &Hedges[40], &Hedges[39]);
	Hedges[39].setAllPtr(&tableVertex[5], &faces[16], &Hedges[49], &Hedges[52], &Hedges[38]);
	Hedges[40].setAllPtr(&tableVertex[5], &faces[6], &Hedges[38], &Hedges[37], &Hedges[41]);
	Hedges[41].setAllPtr(&tableVertex[4], &faces[5], &Hedges[42], &Hedges[44], &Hedges[40]);
	Hedges[42].setAllPtr(&tableVertex[5], &faces[5], &Hedges[41], &Hedges[44], &Hedges[43]);
	Hedges[43].setAllPtr(&tableVertex[9], &faces[14], &Hedges[19], &Hedges[46], &Hedges[42]);
	Hedges[44].setAllPtr(&tableVertex[9], &faces[5], &Hedges[42], &Hedges[41], &Hedges[45]);
	Hedges[45].setAllPtr(&tableVertex[4], &faces[4], &Hedges[5], &Hedges[15], &Hedges[44]);
	Hedges[46].setAllPtr(&tableVertex[5], &faces[14], &Hedges[43], &Hedges[19], &Hedges[47]);
	Hedges[47].setAllPtr(&tableVertex[3], &faces[15], &Hedges[50], &Hedges[48], &Hedges[46]);
	Hedges[48].setAllPtr(&tableVertex[5], &faces[15], &Hedges[47], &Hedges[50], &Hedges[49]);
	Hedges[49].setAllPtr(&tableVertex[2], &faces[16], &Hedges[52], &Hedges[39], &Hedges[44]);
	Hedges[50].setAllPtr(&tableVertex[2], &faces[15], &Hedges[48], &Hedges[47], &Hedges[51]);
	Hedges[51].setAllPtr(&tableVertex[3], &faces[19], &Hedges[58], &Hedges[57], &Hedges[50]);
	Hedges[52].setAllPtr(&tableVertex[8], &faces[16], &Hedges[39], &Hedges[49], &Hedges[53]);
	Hedges[53].setAllPtr(&tableVertex[2], &faces[17], &Hedges[54], &Hedges[33], &Hedges[52]);
	Hedges[54].setAllPtr(&tableVertex[10], &faces[17], &Hedges[33], &Hedges[53], &Hedges[55]);
	Hedges[55].setAllPtr(&tableVertex[2], &faces[18], &Hedges[56], &Hedges[27], &Hedges[54]);
	Hedges[56].setAllPtr(&tableVertex[7], &faces[18], &Hedges[27], &Hedges[55], &Hedges[57]);
	Hedges[57].setAllPtr(&tableVertex[2], &faces[19], &Hedges[51], &Hedges[58], &Hedges[56]);
	Hedges[58].setAllPtr(&tableVertex[7], &faces[19], &Hedges[57], &Hedges[51], &Hedges[59]);
	Hedges[59].setAllPtr(&tableVertex[3], &faces[12], &Hedges[21], &Hedges[22], &Hedges[58]);

	



}

void ViewerWidget::fillFaces()
{

	faces[0].setPointer(&Hedges[0]);
	faces[1].setPointer(&Hedges[1]);
	faces[2].setPointer(&Hedges[7]);
	faces[3].setPointer(&Hedges[11]);
	faces[4].setPointer(&Hedges[5]);
	faces[5].setPointer(&Hedges[41]);
	faces[6].setPointer(&Hedges[37]);
	faces[7].setPointer(&Hedges[3]);
	faces[8].setPointer(&Hedges[31]);
	faces[9].setPointer(&Hedges[9]);
	faces[10].setPointer(&Hedges[25]);
	faces[11].setPointer(&Hedges[13]);
	faces[12].setPointer(&Hedges[21]);
	faces[13].setPointer(&Hedges[17]);
	faces[14].setPointer(&Hedges[19]);
	faces[15].setPointer(&Hedges[47]);
	faces[16].setPointer(&Hedges[39]);
	faces[17].setPointer(&Hedges[33]);
	faces[18].setPointer(&Hedges[27]);
	faces[19].setPointer(&Hedges[51]);

}


//Slots
void ViewerWidget::paintEvent(QPaintEvent* event)
{
	QPainter painter(this);
	QRect area = event->rect();
	painter.drawImage(area, *img, area);
}