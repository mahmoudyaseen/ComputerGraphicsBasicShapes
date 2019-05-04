#if defined(UNICODE) && !defined(_UNICODE)
    #define _UNICODE
#elif defined(_UNICODE) && !defined(UNICODE)
    #define UNICODE
#endif
#include <fstream>
#include <tchar.h>
#include <windows.h>
#include <math.h>
using namespace std;

////////////////////////
void AddMenus(HWND);
HMENU hmenu;
////////////////////////
static fstream  temp ;
////////////////////////

/*  Declare Windows procedure  */
LRESULT CALLBACK WindowProcedure (HWND, UINT, WPARAM, LPARAM);

/*  Make the class name into a global variable  */
TCHAR szClassName[ ] = _T("CodeBlocksWindowsApp");

int WINAPI WinMain (HINSTANCE hThisInstance,
                     HINSTANCE hPrevInstance,
                     LPSTR lpszArgument,
                     int nCmdShow)
{
    temp.open ("temp.txt", ios::out );  /////////////////////////////open temp file
    HWND hwnd;               /* This is the handle for our window */
    MSG messages;            /* Here messages to the application are saved */
    WNDCLASSEX wincl;        /* Data structure for the windowclass */

    /* The Window structure */
    wincl.hInstance = hThisInstance;
    wincl.lpszClassName = szClassName;
    wincl.lpfnWndProc = WindowProcedure;      /* This function is called by windows */
    wincl.style = CS_DBLCLKS;                 /* Catch double-clicks */
    wincl.cbSize = sizeof (WNDCLASSEX);

    /* Use default icon and mouse-pointer */
    wincl.hIcon = LoadIcon (NULL, IDI_APPLICATION);
    wincl.hIconSm = LoadIcon (NULL, IDI_APPLICATION);
    wincl.hCursor = LoadCursor (NULL, IDC_ARROW);
    wincl.lpszMenuName = NULL;                 /* No menu */
    wincl.cbClsExtra = 0;                      /* No extra bytes after the window class */
    wincl.cbWndExtra = 0;                      /* structure or the window instance */
    /* Use Windows's default colour as the background of the window */
    wincl.hbrBackground = (HBRUSH) COLOR_BACKGROUND;

    /* Register the window class, and if it fails quit the program */
    if (!RegisterClassEx (&wincl))
        return 0;

    /* The class is registered, let's create the program*/
    hwnd = CreateWindowEx (
           0,                   /* Extended possibilites for variation */
           szClassName,         /* Classname */
           _T("Code::Blocks Template Windows App"),       /* Title Text */
           WS_OVERLAPPEDWINDOW, /* default window */
           CW_USEDEFAULT,       /* Windows decides the position */
           CW_USEDEFAULT,       /* where the window ends up on the screen */
           544,                 /* The programs width */
           375,                 /* and height in pixels */
           HWND_DESKTOP,        /* The window is a child-window to desktop */
           NULL,                /* No menu */
           hThisInstance,       /* Program Instance handler */
           NULL                 /* No Window Creation data */
           );

    /* Make the window visible on the screen */
    ShowWindow (hwnd, nCmdShow);

    /* Run the message loop. It will run until GetMessage() returns 0 */
    while (GetMessage (&messages, NULL, 0, 0))
    {
        /* Translate virtual-key messages into character messages */
        TranslateMessage(&messages);
        /* Send message to WindowProcedure */
        DispatchMessage(&messages);
    }
    /* The program return-value is 0 - The value that PostQuitMessage() gave */
    temp.close(); ///////////close temp
    return messages.wParam;
}

int Round(double num) { return (int)(num + 0.5); }
void swap(int &x1, int &x2)
{
	int temp = x1;
	x1 = x2;
	x2 = temp;

}
void DrawLineParametric(HDC hdc, int XStart, int YStart, int XEnd, int YEnd, COLORREF R)
{
	int XTemp, YTemp;
	int dx = XEnd - XStart, dy = YEnd - YStart;
	double dt = 1.0 / ( max( abs(dx), abs(dy) ) + 1);
	//dt = 0.0001
	for (double t = 0.0; t < 1; t += dt)
	{
		XTemp = Round(XStart + t * dx);
		YTemp = Round(YStart + t * dy);
		SetPixel(hdc, XTemp, YTemp, R);
	}
}

void DrawLineDDA(HDC hdc, int XStart, int YStart, int XEnd, int YEnd, COLORREF R)
{
	int dx = abs(XEnd - XStart), dy = abs(YEnd - YStart);
	if (dx > dy)
	{
		if (XStart > XEnd)
		{
			swap(XStart, XEnd);
			swap(YStart, YEnd);
		}
		int e = 2 * dy - dx, e1 = 2 * dx, e2 = 2 * dy;
		int XTemp = XStart, YTemp = YStart;
		while (XTemp <= XEnd)
		{
			SetPixel(hdc, XTemp, YTemp, R);
			XTemp++;
			if (e >= 0)
			{
				if (YStart < YEnd)
					YTemp++;
				else
					YTemp--;
				e -= e1;
			}
			e += e2;
		}
	}
	else
	{
		if (YStart > YEnd)
		{
			swap(XStart, XEnd);
			swap(YStart, YEnd);
		}
		int e = 2 * dx - dy, e1 = 2 * dy, e2 = 2 * dx;
		int XTemp = XStart, YTemp = YStart;
		while (YTemp <= YEnd)
		{
			SetPixel(hdc, XTemp, YTemp, R);
			YTemp++;
			if (e >= 0)
			{
				if (XStart < XEnd)
					XTemp++;
				else
					XTemp--;
				e -= e1;
			}
			e += e2;
		}
	}
}

void DrawLineMidPoind(HDC hdc, int XStart, int YStart, int XEnd, int YEnd, COLORREF R)
{
	int dx = XEnd - XStart, dy = YEnd - YStart;
	if (abs(dx) > abs(dy))
	{
		if (XStart > XEnd)
		{
			swap(XStart, XEnd);
			swap(YStart, YEnd);
			dx = -1 * dx;
			dy = -1 * dy;
		}
		if (YStart < YEnd)
		{//slope >= 0 && <=45 || slope += 180 after swap
			int DInitial = dx - 2 * dy;
			int D1 = 2 * dx - 2 * dy; //x++ y--  d<=0 //because we see if we increase y or not when point under or above line
			int D2 = -2 * dy; // x++ d>0
			int XTemp = XStart, YTemp = YStart;
			while (XTemp <= XEnd)
			{
				SetPixel(hdc, XTemp, YTemp, R);
				XTemp++;
				if (DInitial > 0)
				{
					DInitial += D2;
				}
				else
				{
					YTemp++;
					DInitial += D1;
				}
			}
		}
		else
		{	//slope >= 135 && <=180 || slope += 180 after swap
			int DInitial = -1 * dx - 2 * dy;
			int D1 = -2 * dx - 2 * dy; //x++ y--  d<=0 //because we see if we increase y or not when point under or above line
			int D2 = -2 * dy; // x++ d>0
			int XTemp = XStart, YTemp = YStart;
			while (XTemp <= XEnd)
			{
				SetPixel(hdc, XTemp, YTemp, R);
				XTemp++;
				if (DInitial > 0)
				{
					YTemp--;
					DInitial += D1;
				}
				else
				{

					DInitial += D2;
				}
			}

		}
	}
	else
	{
		if (YStart > YEnd)
		{
			swap(XStart, XEnd);
			swap(YStart, YEnd);
			dx = -1 * dx;
			dy = -1 * dy;
		}
		if (XStart < XEnd)
		{//slope >= 45 && <=90 || slope += 180 after swap
			int DInitial = 2 * dx - dy;
			int D1 = 2 * dx - 2 * dy; //x++  d<=0  //swap because we see if we increase x or not when point under or above line
			int D2 = 2 * dx; // y++ x++ d>0
			int XTemp = XStart, YTemp = YStart;
			while (YTemp <= YEnd)
			{
				SetPixel(hdc, XTemp, YTemp, R);
				YTemp++;
				if (DInitial <= 0)
				{
					DInitial += D2;
				}
				else
				{
					XTemp++;
					DInitial += D1;
				}
			}
		}
		else
		{//slope >= 90 && <=135 || slope += 180 after swap
			int DInitial = 2 * dx + dy;
			int D1 = 2 * dx ; //y++  d<=0  //swap because we see if we increase x or not when point under or above line
			int D2 = 2 * dx + 2 * dy; // y++ x-- d>0
			int XTemp = XStart, YTemp = YStart;
			while (YTemp <= YEnd)
			{
				SetPixel(hdc, XTemp, YTemp, R);
				YTemp++;
				if (DInitial > 0)
				{
					DInitial += D1;
				}
				else
				{
					XTemp--;
					DInitial += D2;
				}
			}

		}
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void Draw4Points(HDC hdc, int XC, int YC, int X , int Y , COLORREF R)
{
	SetPixel(hdc, XC + X, YC + Y, R);
	SetPixel(hdc, XC - X, YC + Y, R);
	SetPixel(hdc, XC + X, YC - Y, R);
	SetPixel(hdc, XC - X, YC - Y, R);
}

void swap(int &x1, int &y1, int &x2, int &y2)
{
	swap( x1,  x2);
	swap( y1,  y2);
}

void DrawEllipsDirect(HDC hdc, int XC, int YC, int APointX, int APointY, int BPointX,int  BPointY , COLORREF R)
{
	if (abs(APointX - XC) < abs(APointY - YC) && abs(BPointY - YC) < abs(BPointX - XC))
		swap(APointX, APointY, BPointX, BPointY);
	int A = abs(APointX - XC), B = abs(BPointY - YC);
	int ASqr = A * A, BSqr = B * B;
	double X = A; double Y = 0;
	while (X * BSqr >= Y * ASqr)
	{
		Draw4Points(hdc,  XC,  YC, Round(X), Y,  R);
		Y++;
		X = sqrt((1 - (double)(Y * Y) / BSqr) * ASqr);
	}
	 X = 0;   Y = B;
	while (X * BSqr <= Y * ASqr)
	{
		Draw4Points(hdc, XC, YC, X , Round(Y), R);
		X++;
		Y = sqrt((1 - (double)(X * X) / ASqr) * BSqr);
	}
}

void DrawEllipsParametric(HDC hdc, int XC, int YC, int APointX, int APointY, int BPointX, int  BPointY, COLORREF R)
{
	if (abs(APointX - XC) < abs(APointY - YC) && abs(BPointY - YC) < abs(BPointX - XC))
		swap(APointX, APointY, BPointX, BPointY);
	double x, y;
	int A = abs(APointX - XC), B = abs(BPointY - YC);
	for (double theta = 0; theta < 1.57; theta+=0.001)
	{
		x = Round(A * cos(theta));
		y = Round(B * sin(theta));
		Draw4Points(hdc, XC, YC, x, y, R);
	}
}

void DrawEllipsPolar(HDC hdc, int XC, int YC, int APointX, int APointY, int BPointX, int  BPointY, COLORREF R)
{
	if (abs(APointX - XC) < abs(APointY - YC) && abs(BPointY - YC) < abs(BPointX - XC))
		swap(APointX, APointY, BPointX, BPointY);
	int A = abs(APointX - XC), B = abs(BPointY - YC);
	double x = A, y = 0;
	double theta = 1 / (0.7 * (A + B));
	double costheta = cos(theta), sintheta = sin(theta);
	while (x >= 0)
	{
		double temp = x * costheta - A * y / B * sintheta;
		y = B * x / A * sintheta + y * costheta;
		x = temp;
		Draw4Points(hdc, XC, YC, Round(x), Round(y), R);
	}
}

//if he swap a b

void DrawEllipsMidPoint(HDC hdc, int XC, int YC, int APointX, int APointY, int BPointX, int  BPointY, COLORREF R)
{
	if (abs(APointX - XC) < abs(APointY - YC) && abs(BPointY - YC) < abs(BPointX - XC))
		swap(APointX, APointY, BPointX, BPointY);
	int A = abs(APointX - XC), B = abs(BPointY - YC);
	int ASqr = A * A, BSqr = B * B;
	int dinitial = 4 * BSqr - 4 * ASqr * B +  ASqr, d1 = 12 * BSqr, d2 = 12 * BSqr + 8 * ASqr - 8 * ASqr * B;
	int temp1 = 8 * BSqr, temp2 = 8 * BSqr + 8 * ASqr;
	int X = 0; int Y = B;
	while (X * BSqr <= Y * ASqr)
	{
		if (dinitial <= 0)
		{
			dinitial += d1;
			d1 += temp1;
			d2 += temp1;
			X++;
		}
		else
		{
			dinitial += d2;
			d1 += temp1;
			d2 += temp2;
			X++;
			Y--;
		}
		Draw4Points(hdc, XC, YC, X, Y, R);
	}
	dinitial = 4 * ASqr - 4 * BSqr * A + BSqr; d1 = 12 * ASqr; d2 = 12 * ASqr + 8 * BSqr - 8 * BSqr * A;
	temp1 = 8 * ASqr; temp2 = 8 * BSqr + 8 * ASqr;
	X = A; Y = 0;
	while (X * BSqr >= Y * ASqr)
	{
		if (dinitial <= 0)
		{
			dinitial += d1;
			d1 += temp1;
			d2 += temp1;
			Y++;
		}
		else
		{
			dinitial += d2;
			d1 += temp1;
			d2 += temp2;
			X--;
			Y++;
		}
		Draw4Points(hdc, XC, YC, X, Y, R);
	}

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct Pair
{
	int first;
	int second;
};
void DrawRectangleWindow(HDC hdc, int XStart, int YStart, int XEnd, int YEnd, COLORREF R)
{
	DrawLineMidPoind(hdc, XStart, YStart, XEnd, YStart, R);
	DrawLineMidPoind(hdc, XStart, YStart, XStart, YEnd, R);
	DrawLineMidPoind(hdc, XEnd, YStart, XEnd, YEnd, R);
	DrawLineMidPoind(hdc, XStart, YEnd, XEnd, YEnd, R);
}
bool IfInLeftLine(int XStart, Pair point) { return point.first > XStart; }
bool IfInRightLine(int XEnd, Pair point) { return point.first < XEnd; }
bool IfInUpLine(int YStart, Pair point) { return point.second > YStart; }
bool IfInDownLine(int YEnd, Pair point) { return point.second < YEnd; }
Pair GetIntersectVertical(Pair point1, Pair point2, int X)
{
	Pair intersect;
	intersect.first = X;
	intersect.second = Round(point1.second + (double)(X - point1.first) * (point2.second - point1.second) / (point2.first - point1.first));
	return intersect;
}
Pair GetIntersectHorisontal(Pair point1, Pair point2, int Y)
{
	Pair intersect;
	intersect.second = Y;
	intersect.first = Round(point1.first + (double)(Y - point1.second) * (point2.first - point1.first) / (point2.second - point1.second));
	return intersect;
}
Pair* VerticalClipping(int X, Pair * points , int size , int type , int &OutSize)
{//1 for left 2 for right
	Pair* NewPoints = new Pair[200];
	OutSize = 0;
	for (int i = 0; i < size; i++)
	{
		bool point1 ;
		bool point2 ;
		if (type == 1)
		{
			point1 = IfInLeftLine(X, points[i]);
			point2 = IfInLeftLine(X, points[(i + 1) % size]);
		}
		else
		{
			point1 = IfInRightLine(X, points[i]);
			point2 = IfInRightLine(X, points[(i + 1) % size]);
		}
		if (point1 && !point2)
		{
			NewPoints[OutSize++] = GetIntersectVertical(points[i], points[(i + 1) % size], X);
		}
		else if (!point1 && point2)
		{
			NewPoints[OutSize++] = GetIntersectVertical(points[i], points[(i + 1) % size], X);
			NewPoints[OutSize++] = (points[(i + 1) % size]);
		}
		else if (point1 && point2)
		{
			NewPoints[OutSize++] = (points[(i + 1) % size]);
		}
	}
	return NewPoints;
}
Pair* HorizontalClipping(int Y, Pair * points, int size, int type, int &OutSize)
{//1 for up 2 for down
	Pair* NewPoints = new Pair[200];
	OutSize = 0;
	for (int i = 0; i < size; i++)
	{
		bool point1;
		bool point2;
		if (type == 1)
		{
			point1 = IfInUpLine(Y, points[i]);
			point2 = IfInUpLine(Y, points[(i + 1) % size]);
		}
		else
		{
			point1 = IfInDownLine(Y, points[i]);
			point2 = IfInDownLine(Y, points[(i + 1) % size]);
		}
		if (point1 && !point2)
		{
			NewPoints[OutSize++] = GetIntersectHorisontal(points[i], points[(i + 1) % size], Y);

		}
		else if (!point1 && point2)
		{
			NewPoints[OutSize++] = GetIntersectHorisontal(points[i], points[(i + 1) % size], Y);
			NewPoints[OutSize++] = (points[(i + 1) % size]);
		}
		else if (point1 && point2)
		{
			NewPoints[OutSize++] = (points[(i + 1) % size]);
		}
	}
	return NewPoints;
}

void PolygonClipping(HDC hdc, int XStart, int YStart, int XEnd, int YEnd, Pair* points , int size, COLORREF R)
{
    if(XStart > XEnd) swap(XStart , XEnd);
    if(YStart > YEnd) swap(YStart , YEnd);
	for (int i = 0; i < size; i++)
	{
		int xs = points[i].first, ys = points[i].second, xe = points[(i + 1) % size].first, ye = points[(i + 1) % size].second;
		DrawLineMidPoind(hdc, xs, ys, xe, ye, RGB(0,0,0));
	}
	DrawRectangleWindow(hdc, XStart, YStart, XEnd, YEnd, RGB(255, 0, 0));
	int newsize;
	Pair* newpoints;
	newpoints = VerticalClipping(XStart, points, size, 1, newsize);//left
	newpoints = VerticalClipping(XEnd, newpoints, newsize, 2, newsize);//right
	newpoints = HorizontalClipping(YStart, newpoints, newsize, 1, newsize);//up
	newpoints = HorizontalClipping(YEnd, newpoints, newsize, 2, newsize);//down
	for (int i = 0; i < newsize; i++)
	{
		int xs = newpoints[i].first, ys = newpoints[i].second, xe = newpoints[(i+1)%newsize].first, ye = newpoints[(i + 1) % newsize].second;
		DrawLineMidPoind(hdc, xs, ys, xe, ye, R);
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int multi(double *TVector, int **Matrix, int* XVector , int size)
{
	int *res = new int[size];
	for (int i = 0; i < size; i++)
		res[i] = 0;
	for (int i = 0; i < size; i++)
		for (int j = 0; j < size; j++)
			res[i] += Matrix[i][j] * XVector[j];
	double Total = 0;
	for (int i = 0; i < size; i++)
		Total += res[i] * TVector[i];
	return Total;
}

void PassSpline(HDC hdc, int XStart, int YStart, int XMid, int YMid, int XEnd, int YEnd, COLORREF R)
{
	double *TVector = new double[3];
	int **Matrix = new int*[3];
	for (int i = 0; i < 3; i++)
		Matrix[i] = new int[3];
	int *XVector = new int[3];
	int *YVector = new int[3];
	Matrix[0][0] = 1; Matrix[0][1] = 0; Matrix[0][2] = 0;
	Matrix[1][0] = -3; Matrix[1][1] = 4; Matrix[1][2] = -1;
	Matrix[2][0] = 2; Matrix[2][1] = -4; Matrix[2][2] = 2;
	XVector[0] = XStart; XVector[1] = XMid; XVector[2] = XEnd;
	YVector[0] = YStart; YVector[1] = YMid; YVector[2] = YEnd;
	for (double t = 0.0; t < 1; t+=0.0001)
	{
		TVector[0] = 1; TVector[1] = t; TVector[2] = t * t;
		int x = Round(multi(TVector, Matrix, XVector, 3));
		int y = Round(multi(TVector, Matrix, YVector, 3));
		SetPixel(hdc, x, y, R);
	}
}

void GideSpline(HDC hdc, int XStart, int YStart, int XMid, int YMid, int XEnd, int YEnd, COLORREF R)
{
	double *TVector = new double[3];
	int **Matrix = new int*[3];
	for (int i = 0; i < 3; i++)
		Matrix[i] = new int[3];
	int *XVector = new int[3];
	int *YVector = new int[3];
	Matrix[0][0] = 1; Matrix[0][1] = 0; Matrix[0][2] = 0;
	Matrix[1][0] = -2; Matrix[1][1] = 2; Matrix[1][2] = 0;
	Matrix[2][0] = 1; Matrix[2][1] = -2; Matrix[2][2] = 1;
	XVector[0] = XStart; XVector[1] = XMid; XVector[2] = XEnd;
	YVector[0] = YStart; YVector[1] = YMid; YVector[2] = YEnd;
	for (double t = 0.0; t < 1; t += 0.0001)
	{
		TVector[0] = 1; TVector[1] = t; TVector[2] = t * t;
		int x = Round(multi(TVector, Matrix, XVector, 3));
		int y = Round(multi(TVector, Matrix, YVector, 3));
		SetPixel(hdc, x, y, R);
	}
}

void HermiteCurve(HDC hdc, int XStart, int YStart, int UStart, int VStart, int XEnd, int YEnd, int UEnd, int VEnd, COLORREF R)
{
	double *TVector = new double[4];
	int **Matrix = new int*[4];
	for (int i = 0; i < 4; i++)
		Matrix[i] = new int[4];
	int *XVector = new int[4];
	int *YVector = new int[4];
	Matrix[0][0] = 1; Matrix[0][1] = 0; Matrix[0][2] = 0; Matrix[0][3] = 0;
	Matrix[1][0] = 0; Matrix[1][1] = 1; Matrix[1][2] = 0; Matrix[1][3] = 0;
	Matrix[2][0] = -3; Matrix[2][1] = -2; Matrix[2][2] = 3; Matrix[2][3] = -1;
	Matrix[3][0] = 2; Matrix[3][1] = 1; Matrix[3][2] = -2; Matrix[3][3] = 1;
	XVector[0] = XStart; XVector[1] = UStart; XVector[2] = XEnd; XVector[3] = UEnd;
	YVector[0] = YStart; YVector[1] = VStart; YVector[2] = YEnd; YVector[3] = VEnd;
	for (double t = 0.0; t < 1; t += 0.0001)
	{
		TVector[0] = 1; TVector[1] = t; TVector[2] = t * t; TVector[3] = t * t * t;
		int x = Round(multi(TVector, Matrix, XVector, 4));
		int y = Round(multi(TVector, Matrix, YVector, 4));
		SetPixel(hdc, x, y, R);
	}

}

void BeizerCurve(HDC hdc, int X1, int Y1, int X2, int Y2, int X3, int Y3, int X4, int Y4, COLORREF R)
{
	double *TVector = new double[4];
	int **Matrix = new int*[4];
	for (int i = 0; i < 4; i++)
		Matrix[i] = new int[4];
	int *XVector = new int[4];
	int *YVector = new int[4];
	Matrix[0][0] = 1; Matrix[0][1] = 0; Matrix[0][2] = 0; Matrix[0][3] = 0;
	Matrix[1][0] = -3; Matrix[1][1] = 3; Matrix[1][2] = 0; Matrix[1][3] = 0;
	Matrix[2][0] = 3; Matrix[2][1] = -6; Matrix[2][2] = 3; Matrix[2][3] = 0;
	Matrix[3][0] = -1; Matrix[3][1] = 3; Matrix[3][2] = -3; Matrix[3][3] = 1;
	XVector[0] = X1; XVector[1] = X2; XVector[2] = X3; XVector[3] = X4;
	YVector[0] = Y1; YVector[1] = Y2; YVector[2] = Y3; YVector[3] = Y4;
	for (double t = 0.0; t < 1; t += 0.0001)
	{
		TVector[0] = 1; TVector[1] = t; TVector[2] = t * t; TVector[3] = t * t * t;
		int x = Round(multi(TVector, Matrix, XVector, 4));
		int y = Round(multi(TVector, Matrix, YVector, 4));
		SetPixel(hdc, x, y, R);
	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void load(HDC hdc)
{
    temp.close();
    remove("temp.txt");
    CopyFile("seen.txt", "temp.txt", TRUE);
    temp.open ("temp.txt", ios::in );
    int XStart, YStart,  XMid,  YMid, XEnd, YEnd,  UStart,  VStart,   UEnd,  VEnd, Counter = 0 , count = 0;
    static Pair* points = new Pair[200];
	Pair point;
    char choice[200];
    while(!temp.eof())
    {
        temp.getline(choice,200,' ');
        if(strcmp(choice , "ParametricLine") == 0 )
        {
            temp.read((char*)&XStart,sizeof(XStart));
            temp.read((char*)&YStart,sizeof(YStart));
            temp.read((char*)&XEnd,sizeof(XEnd));
            temp.read((char*)&YEnd,sizeof(YEnd));
            DrawLineParametric(hdc, XStart, YStart, XEnd, YEnd, RGB(255, 0, 0));
        }
        else if(strcmp(choice , "DDALine") == 0 )
        {
            temp.read((char*)& XStart,sizeof(XStart));
            temp.read((char*)& YStart,sizeof(YStart));
            temp.read((char*)& XEnd,sizeof(XEnd));
            temp.read((char*)& YEnd,sizeof(YEnd));
            DrawLineDDA(hdc, XStart, YStart, XEnd, YEnd, RGB(255, 0, 0));
        }
        else if(strcmp(choice , "MidPointLine") == 0)
        {
            temp.read((char*)& XStart,sizeof(XStart));
            temp.read((char*)& YStart,sizeof(YStart));
            temp.read((char*)& XEnd,sizeof(XEnd));
            temp.read((char*)& YEnd,sizeof(YEnd));
            DrawLineMidPoind(hdc, XStart, YStart, XEnd, YEnd, RGB(255, 0, 0));
        }
        else if(strcmp(choice , "DirectEllipse") == 0)
        {
            temp.read((char*)& XStart,sizeof(XStart));
            temp.read((char*)& YStart,sizeof(YStart));
            temp.read((char*)& UStart,sizeof(UStart));
            temp.read((char*)& VStart,sizeof(VStart));
            temp.read((char*)& XEnd,sizeof(XEnd));
            temp.read((char*)& YEnd,sizeof(YEnd));
            DrawEllipsDirect(hdc, XStart, YStart, UStart, VStart, XEnd, YEnd, RGB(255, 0, 0));
        }
        else if(strcmp(choice , "PolarEllipse") == 0)
        {
            temp.read((char*)& XStart,sizeof(XStart));
            temp.read((char*)& YStart,sizeof(YStart));
            temp.read((char*)& UStart,sizeof(UStart));
            temp.read((char*)& VStart,sizeof(VStart));
            temp.read((char*)& XEnd,sizeof(XEnd));
            temp.read((char*)& YEnd,sizeof(YEnd));
            DrawEllipsPolar(hdc, XStart, YStart, UStart, VStart, XEnd, YEnd, RGB(255, 0, 0));
        }
        else if(strcmp(choice , "ParametricEllipse") == 0)
        {
            temp.read((char*)& XStart,sizeof(XStart));
            temp.read((char*)& YStart,sizeof(YStart));
            temp.read((char*)& UStart,sizeof(UStart));
            temp.read((char*)& VStart,sizeof(VStart));
            temp.read((char*)& XEnd,sizeof(XEnd));
            temp.read((char*)& YEnd,sizeof(YEnd));
            DrawEllipsParametric(hdc, XStart, YStart, UStart, VStart, XEnd, YEnd, RGB(255, 0, 0));
        }
        else if(strcmp(choice , "MidPointEllipse") == 0)
        {
            temp.read((char*)& XStart,sizeof(XStart));
            temp.read((char*)& YStart,sizeof(YStart));
            temp.read((char*)& UStart,sizeof(UStart));
            temp.read((char*)& VStart,sizeof(VStart));
            temp.read((char*)& XEnd,sizeof(XEnd));
            temp.read((char*)& YEnd,sizeof(YEnd));
            DrawEllipsMidPoint(hdc, XStart, YStart, UStart, VStart, XEnd, YEnd, RGB(255, 0, 0));
        }
        else if(strcmp(choice , "splinesPass") == 0)
        {
            temp.read((char*)& XStart,sizeof(XStart));
            temp.read((char*)& YStart,sizeof(YStart));
            temp.read((char*)& XMid,sizeof(XMid));
            temp.read((char*)& YMid,sizeof(YMid));
            temp.read((char*)& XEnd,sizeof(XEnd));
            temp.read((char*)& YEnd,sizeof(YEnd));
            PassSpline(hdc, XStart, YStart, XMid, YMid, XEnd, YEnd, RGB(255, 0, 0));
        }
        else if(strcmp(choice , "splinesGide") == 0)
        {
            temp.read((char*)& XStart,sizeof(XStart));
            temp.read((char*)& YStart,sizeof(YStart));
            temp.read((char*)& XMid,sizeof(XMid));
            temp.read((char*)& YMid,sizeof(YMid));
            temp.read((char*)& XEnd,sizeof(XEnd));
            temp.read((char*)& YEnd,sizeof(YEnd));
            GideSpline(hdc, XStart, YStart, XMid, YMid, XEnd, YEnd, RGB(255, 0, 0));
        }
        else if(strcmp(choice , "bezier") == 0)
        {
            temp.read((char*)& XStart,sizeof(XStart));
            temp.read((char*)& YStart,sizeof(YStart));
            temp.read((char*)& UStart,sizeof(UStart));
            temp.read((char*)& VStart,sizeof(VStart));
            temp.read((char*)& XEnd,sizeof(XEnd));
            temp.read((char*)& YEnd,sizeof(YEnd));
            temp.read((char*)& UEnd,sizeof(UEnd));
            temp.read((char*)& VEnd,sizeof(VEnd));
            BeizerCurve(hdc, XStart, YStart, UStart, VStart, XEnd, YEnd, UEnd, VEnd, RGB(255, 0, 0));
        }
        else if(strcmp(choice , "hermite") == 0)
        {
            temp.read((char*)& XStart,sizeof(XStart));
            temp.read((char*)& YStart,sizeof(YStart));
            temp.read((char*)& UStart,sizeof(UStart));
            temp.read((char*)& VStart,sizeof(VStart));
            temp.read((char*)& XEnd,sizeof(XEnd));
            temp.read((char*)& YEnd,sizeof(YEnd));
            temp.read((char*)& UEnd,sizeof(UEnd));
            temp.read((char*)& VEnd,sizeof(VEnd));
            HermiteCurve(hdc, XStart, YStart, UStart, VStart, XEnd, YEnd, UEnd, VEnd, RGB(255, 0, 0));
        }
        else if(strcmp(choice , "clipping") == 0)
        {
            temp.read((char*)& XStart,sizeof(XStart));
            temp.read((char*)& YStart,sizeof(YStart));
            temp.read((char*)& XEnd,sizeof(XEnd));
            temp.read((char*)& YEnd,sizeof(YEnd));
            temp.read((char*)& Counter,sizeof(Counter));
            for(int i = 0 ; i < Counter ; i++)
            {
                temp.read((char*)& points[i].first,sizeof(points[i].first));
                temp.read((char*)& points[i].second,sizeof(points[i].second));
            }
            PolygonClipping(hdc, XStart, YStart, XEnd, YEnd, points, Counter, RGB(255, 0, 0));
        }

    }
    temp.close();
    temp.open ("temp.txt", ios::out | ios::app );
}

void save()
{
    temp.close();
    remove("seen.txt");
    CopyFile("temp.txt", "seen.txt", TRUE);
    temp.open ("temp.txt", ios::out | ios::app );
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*  This function is called by the Windows function DispatchMessage()  */

LRESULT CALLBACK WindowProcedure (HWND hwnd, UINT message, WPARAM wParam, LPARAM lParam)
{
	static int XStart, YStart,  XMid,  YMid, XEnd, YEnd,  UStart,  VStart,   UEnd,  VEnd, Counter = 0 , count = 0;
	HDC hdc = GetDC(hwnd);
	static Pair* points = new Pair[200];
	Pair point;
	static string choice;

    switch (message)                  /* handle the messages */
    {
        case WM_COMMAND:
            {
                switch(wParam)
                {
                case 1:
                    choice = "ParametricLine";
                    break;
                case 2:
                    choice = "DDALine";
                    break;
                case 3:
                    choice = "MidPointLine";
                    break;
                case 4:
                    choice = "DirectEllipse";
                    break;
                case 5:
                    choice = "PolarEllipse";
                    break;
                case 6:
                    choice = "ParametricEllipse";
                    break;
                case 7:
                    choice = "MidPointEllipse";
                    break;
                case 8:
                    choice = "splinesPass";
                    break;
                case 9:
                    choice = "splinesGide";
                    break;
                case 10:
                    choice = "bezier";
                    break;
                case 11:
                    choice = "hermite";
                    break;
                case 12:
                    choice = "Save";
                    save();
                    break;
                case 13:
                    choice = "Load";
                    load(hdc);
                    break;
                case 14:
                    choice = "clipping";
                    break;
                }
            }
            break;
        case WM_LBUTTONUP:
            if(choice == "ParametricLine" || choice == "DDALine" || choice == "MidPointLine")
            {
                if (Counter == 0)
                {
                    XStart = LOWORD(lParam);
                    YStart = HIWORD(lParam);
                    //SetPixel(hdc, XStart, YStart, RGB(255, 0, 0));
                    Counter++;
                }
                else
                {
                    XEnd = LOWORD(lParam);
                    YEnd = HIWORD(lParam);
                    //SetPixel(hdc, XEnd, YEnd, RGB(255, 0, 0));
                    Counter--;
                    if(choice == "ParametricLine")
                    {
                        temp << "ParametricLine ";
                        temp.write((char*)& XStart,sizeof(XStart));
                        temp.write((char*)& YStart,sizeof(YStart));
                        temp.write((char*)& XEnd,sizeof(XEnd));
                        temp.write((char*)& YEnd,sizeof(YEnd));
                        DrawLineParametric(hdc, XStart, YStart, XEnd, YEnd, RGB(255, 0, 0));
                    }
                    else if(choice == "DDALine")
                    {
                        temp << "DDALine ";
                        temp.write((char*)& XStart,sizeof(XStart));
                        temp.write((char*)& YStart,sizeof(YStart));
                        temp.write((char*)& XEnd,sizeof(XEnd));
                        temp.write((char*)& YEnd,sizeof(YEnd));
                        DrawLineDDA(hdc, XStart, YStart, XEnd, YEnd, RGB(255, 0, 0));
                    }
                    else if(choice == "MidPointLine")
                    {
                        temp << "MidPointLine ";
                        temp.write((char*)& XStart,sizeof(XStart));
                        temp.write((char*)& YStart,sizeof(YStart));
                        temp.write((char*)& XEnd,sizeof(XEnd));
                        temp.write((char*)& YEnd,sizeof(YEnd));
                        DrawLineMidPoind(hdc, XStart, YStart, XEnd, YEnd, RGB(255, 0, 0));
                    }
                }
            }
            else if(choice == "DirectEllipse" || choice == "PolarEllipse" || choice == "ParametricEllipse" || choice == "MidPointEllipse")
            {
                if (Counter == 0)
                {
                    XStart = LOWORD(lParam);
                    YStart = HIWORD(lParam);
                    SetPixel(hdc, XStart, YStart, RGB(255, 0, 0));
                    Counter++;
                }
                else if (Counter == 1)
                {
                    UStart = LOWORD(lParam);
                    VStart = HIWORD(lParam);
                    SetPixel(hdc, XEnd, YEnd, RGB(255, 0, 0));
                    Counter++;
                }
                else if (Counter == 2)
                {
                    XEnd = LOWORD(lParam);
                    YEnd = HIWORD(lParam);
                    SetPixel(hdc, XEnd, YEnd, RGB(255, 0, 0));
                    if(choice == "DirectEllipse")
                    {
                        temp << "DirectEllipse ";
                        temp.write((char*)& XStart,sizeof(XStart));
                        temp.write((char*)& YStart,sizeof(YStart));
                        temp.write((char*)& UStart,sizeof(UStart));
                        temp.write((char*)& VStart,sizeof(VStart));
                        temp.write((char*)& XEnd,sizeof(XEnd));
                        temp.write((char*)& YEnd,sizeof(YEnd));
                        DrawEllipsDirect(hdc, XStart, YStart, UStart, VStart, XEnd, YEnd, RGB(255, 0, 0));
                    }
                    else if(choice == "PolarEllipse")
                    {
                        temp << "PolarEllipse ";
                        temp.write((char*)& XStart,sizeof(XStart));
                        temp.write((char*)& YStart,sizeof(YStart));
                        temp.write((char*)& UStart,sizeof(UStart));
                        temp.write((char*)& VStart,sizeof(VStart));
                        temp.write((char*)& XEnd,sizeof(XEnd));
                        temp.write((char*)& YEnd,sizeof(YEnd));
                        DrawEllipsPolar(hdc, XStart, YStart, UStart, VStart, XEnd, YEnd, RGB(255, 0, 0));
                    }
                    else if(choice == "ParametricEllipse")
                    {
                        temp << "ParametricEllipse ";
                        temp.write((char*)& XStart,sizeof(XStart));
                        temp.write((char*)& YStart,sizeof(YStart));
                        temp.write((char*)& UStart,sizeof(UStart));
                        temp.write((char*)& VStart,sizeof(VStart));
                        temp.write((char*)& XEnd,sizeof(XEnd));
                        temp.write((char*)& YEnd,sizeof(YEnd));
                        DrawEllipsParametric(hdc, XStart, YStart, UStart, VStart, XEnd, YEnd, RGB(255, 0, 0));
                    }
                    else if(choice == "MidPointEllipse")
                    {
                        temp << "MidPointEllipse ";
                        temp.write((char*)& XStart,sizeof(XStart));
                        temp.write((char*)& YStart,sizeof(YStart));
                        temp.write((char*)& UStart,sizeof(UStart));
                        temp.write((char*)& VStart,sizeof(VStart));
                        temp.write((char*)& XEnd,sizeof(XEnd));
                        temp.write((char*)& YEnd,sizeof(YEnd));
                        DrawEllipsMidPoint(hdc, XStart, YStart, UStart, VStart, XEnd, YEnd, RGB(255, 0, 0));
                    }
                    Counter = 0;
                }
            }
            else if(choice == "splinesPass" || choice == "splinesGide")
            {
                if (Counter == 0)
                {
                    XStart = LOWORD(lParam);
                    YStart = HIWORD(lParam);
                    //SetPixel(hdc, XStart, YStart, RGB(255, 0, 0));
                    Counter++;
                }
                else if (Counter == 1)
                {
                    XMid = LOWORD(lParam);
                    YMid = HIWORD(lParam);
                    //SetPixel(hdc, XEnd, YEnd, RGB(255, 0, 0));
                    Counter++;
                }
                else if(Counter == 2)
                {
                    XEnd = LOWORD(lParam);
                    YEnd = HIWORD(lParam);
                    //SetPixel(hdc, XEnd, YEnd, RGB(255, 0, 0));
                    if(choice == "splinesPass")
                    {
                        temp << "splinesPass ";
                        temp.write((char*)& XStart,sizeof(XStart));
                        temp.write((char*)& YStart,sizeof(YStart));
                        temp.write((char*)& XMid,sizeof(XMid));
                        temp.write((char*)& YMid,sizeof(YMid));
                        temp.write((char*)& XEnd,sizeof(XEnd));
                        temp.write((char*)& YEnd,sizeof(YEnd));
                        PassSpline(hdc, XStart, YStart, XMid, YMid, XEnd, YEnd, RGB(255, 0, 0));
                    }
                    else if(choice == "splinesGide")
                    {
                        temp << "splinesGide ";
                        temp.write((char*)& XStart,sizeof(XStart));
                        temp.write((char*)& YStart,sizeof(YStart));
                        temp.write((char*)& XMid,sizeof(XMid));
                        temp.write((char*)& YMid,sizeof(YMid));
                        temp.write((char*)& XEnd,sizeof(XEnd));
                        temp.write((char*)& YEnd,sizeof(YEnd));
                        GideSpline(hdc, XStart, YStart, XMid, YMid, XEnd, YEnd, RGB(255, 0, 0));
                    }
                    Counter = 0;
                }
            }
            else if(choice == "bezier" || choice == "hermite")
            {
                if (Counter == 0)
                {
                    XStart = LOWORD(lParam);
                    YStart = HIWORD(lParam);
                    //SetPixel(hdc, XStart, YStart, RGB(255, 0, 0));
                    Counter++;
                }
                else if (Counter == 1)
                {
                    UStart = LOWORD(lParam);
                    VStart = HIWORD(lParam);
                    //SetPixel(hdc, XEnd, YEnd, RGB(255, 0, 0));
                    Counter++;
                }
                else if (Counter == 2)
                {
                    XEnd = LOWORD(lParam);
                    YEnd = HIWORD(lParam);
                    //SetPixel(hdc, XEnd, YEnd, RGB(255, 0, 0));
                    Counter++;
                }
                else if (Counter == 3)
                {
                    UEnd = LOWORD(lParam);
                    VEnd = HIWORD(lParam);
                    //SetPixel(hdc, XEnd, YEnd, RGB(255, 0, 0));
                    if(choice == "bezier")
                    {
                        temp << "bezier ";
                        temp.write((char*)& XStart,sizeof(XStart));
                        temp.write((char*)& YStart,sizeof(YStart));
                        temp.write((char*)& UStart,sizeof(UStart));
                        temp.write((char*)& VStart,sizeof(VStart));
                        temp.write((char*)& XEnd,sizeof(XEnd));
                        temp.write((char*)& YEnd,sizeof(YEnd));
                        temp.write((char*)& UEnd,sizeof(UEnd));
                        temp.write((char*)& VEnd,sizeof(VEnd));
                        BeizerCurve(hdc, XStart, YStart, UStart, VStart, XEnd, YEnd, UEnd, VEnd, RGB(255, 0, 0));
                    }
                    else if(choice == "hermite")
                    {
                        temp << "hermite ";
                        temp.write((char*)& XStart,sizeof(XStart));
                        temp.write((char*)& YStart,sizeof(YStart));
                        temp.write((char*)& UStart,sizeof(UStart));
                        temp.write((char*)& VStart,sizeof(VStart));
                        temp.write((char*)& XEnd,sizeof(XEnd));
                        temp.write((char*)& YEnd,sizeof(YEnd));
                        temp.write((char*)& UEnd,sizeof(UEnd));
                        temp.write((char*)& VEnd,sizeof(VEnd));
                        HermiteCurve(hdc, XStart, YStart, UStart, VStart, XEnd, YEnd, UEnd, VEnd, RGB(255, 0, 0));
                    }
                    Counter = 0;
                }
            }
            else if (choice == "clipping")
            {
                XStart = LOWORD(lParam);
            	YStart = HIWORD(lParam);
            	point.first = XStart; point.second = YStart;
            	points[Counter++] = point;
            }
            break;
        case WM_RBUTTONUP:
            if(choice == "clipping")
            {
                if (count == 0)
                {
                    XStart = LOWORD(lParam);
                    YStart = HIWORD(lParam);
                    count++;
                }
                else
                {
                    XEnd = LOWORD(lParam);
                    YEnd = HIWORD(lParam);
                    temp << "clipping ";
                    temp.write((char*)& XStart,sizeof(XStart));
                    temp.write((char*)& YStart,sizeof(YStart));
                    temp.write((char*)& XEnd,sizeof(XEnd));
                    temp.write((char*)& YEnd,sizeof(YEnd));
                    temp.write((char*)& Counter,sizeof(Counter));
                    for(int i = 0 ; i < Counter ; i++)
                    {
                        temp.write((char*)& points[i].first,sizeof(points[i].first));
                        temp.write((char*)& points[i].second,sizeof(points[i].second));
                    }
                    PolygonClipping(hdc, XStart, YStart, XEnd, YEnd, points, Counter, RGB(255, 0, 0));
                    count = 0;
                    Counter = 0;
                    Pair* points = new Pair[200];
                }
            }
            break;
        case WM_CREATE:
            AddMenus(hwnd);
            break;
        case WM_DESTROY:
            PostQuitMessage (0);       /* send a WM_QUIT to the message queue */
            break;
        default:                      /* for messages that we don't deal with */
            return DefWindowProc (hwnd, message, wParam, lParam);
    }

    return 0;
}

void AddMenus(HWND hwnd)
{
    hmenu = CreateMenu();
    HMENU lineMenu = CreateMenu();
    HMENU ellipse = CreateMenu();
    HMENU clipping = CreateMenu();
    HMENU curves = CreateMenu();
    HMENU store = CreateMenu();

    AppendMenu(lineMenu , MF_STRING , 1 , "Parametric"); //id = 1
    AppendMenu(lineMenu , MF_STRING , 2 , "DDA"); //id = 2
    AppendMenu(lineMenu , MF_STRING , 3 , "MidPoint"); //id = 3

    AppendMenu(ellipse , MF_STRING , 4 , "Direct"); //id = 1
    AppendMenu(ellipse , MF_STRING , 5 , "polar"); //id = 2
    AppendMenu(ellipse , MF_STRING , 6 , "Parametric"); //id = 2
    AppendMenu(ellipse , MF_STRING , 7 , "midpoint"); //id = 3

    AppendMenu(curves , MF_STRING , 8 , "splinesPass"); //id = 1
    AppendMenu(curves , MF_STRING , 9 , "splinesGide"); //id = 1
    AppendMenu(curves , MF_STRING , 10 , "bezier"); //id = 2
    AppendMenu(curves , MF_STRING , 11 , "hermite"); //id = 3

    AppendMenu(store , MF_STRING , 12 , "Save"); //id = 1
    AppendMenu(store , MF_STRING , 13 , "Load"); //id = 1


    AppendMenu(hmenu , MF_POPUP , (UINT_PTR)lineMenu , "Line");
    AppendMenu(hmenu , MF_POPUP , (UINT_PTR)ellipse , "Ellipse");
    AppendMenu(hmenu , MF_STRING , 14 , "clipping"); //id = 1
    AppendMenu(hmenu , MF_POPUP , (UINT_PTR)curves , "Curves");
    AppendMenu(hmenu , MF_POPUP , (UINT_PTR)store , "Store");
    SetMenu(hwnd , hmenu);

}
