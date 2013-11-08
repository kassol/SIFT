#pragma once
#include "sift.h"

struct Pt
{
	Pt(){x = 0; y = 0;}
	Pt(double dx, double dy){x = dx; y = dy;}
	Pt(int nx, int ny){x = (double)nx; y = (double)ny;}
	double x;
	double y;
};

struct RBTree
{
	kd_node* node;
	Keypoint* feature;
	int num;
	RBTree(){node = NULL; feature = NULL;num = 0;}
	RBTree(kd_node* n, Keypoint* feat, Keypoint** nb){node = n; feature = feat; num = 0;}
};

struct Rec
{
	Rec(){left = 0; right = 0; top = 0; bottom = 0;}
	Rec(double l, double r, double t, double b){left = l; right = r; top = t; bottom = b;}
	bool Intersects(Rec& rect)
	{
		if ((left-rect.right)*(right-rect.left) < 0 && (top - rect.bottom)*(bottom - rect.top) < 0)
		{
			return true;
		}
		return false;
	}
	Rec Intersected(Rec& rect)
	{
		if (this->Intersects(rect))
		{
			return Rec(left>rect.left?left:rect.left, right<rect.right?right:rect.right, top>rect.top?top:rect.top, bottom<rect.bottom?bottom:rect.bottom);
		}
		return Rec(0, 0, 0, 0);
	}
	Rec Union(Rec& rect)
	{
		Rec temRec(left<rect.left?left:rect.left, right>rect.right?right:rect.right, top<rect.top?top:rect.top, bottom>rect.bottom?bottom:rect.bottom);
	}
	bool IsEmpty()
	{
		if (right-left == 0 || top-bottom == 0)
		{
			return true;
		}
		return false;
	}
	void extend()
	{
		left = floor(left);
		right = ceil(right);
		top = floor(top);
		bottom = ceil(bottom);
	}
	double left;
	double right;
	double top;
	double bottom;
};


class match
{
public:
	match(char* l, char* r);
	~match(void);
	void domatch(std::vector<SamePoint>& resultData);

private:
	char* m_szPathNameL;
	char* m_szPathNameR;
	Sift sl;
	Sift sr;
	IImage* m_pImage;
	IImage* m_pImage2;
};

